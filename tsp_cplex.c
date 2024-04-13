#include "tsp_cplex.h"
#include "tsp.h"
#include "util.h"
#include <stdio.h>
#include <string.h>

int xpos(int i, int j, struct tsp* tsp)
{
	if (i == j) {
		printf("i == j\n");
		return -1;
	}
	if (i > j)
		return xpos(j, i, tsp);
	return i * tsp->nnodes + j - ((i + 1) * (i + 2)) / 2;
}

int tsp_build_lpmodel(struct tsp* tsp, CPXENVptr env, CPXLPptr lp)
{
	char binary = 'B';
	char* col_name = (char*)calloc(100, sizeof(char));

	// add variables to cplex
	for (int i = 0; i < tsp->nnodes - 1; i++) {
		for (int j = i + 1; j < tsp->nnodes; j++) {
			double cost = tsp->cost_matrix[flatten_coords(i, j, tsp->nnodes)];
			double lb = 0;
			double ub = 1;
			sprintf(col_name, "x(%d,%d)", i + 1, j + 1);
			int err;
			if ((err = CPXnewcols(env, lp, 1, &cost, &lb, &ub, &binary, &col_name))) {
				printf("Error adding variable: %d\n", err);
				return -1;
			}
		}
	}

	free(col_name);

	// add contraints to cplex
	double rhs = 2;
	char sense = 'E';
	int start_row_coefficients = 0; // it would be an array containing the start of each row's coefficients in
					// vars and coeffs
	int* vars = (int*)calloc(tsp->nnodes - 1, sizeof(int));
	double* coeffs = (double*)calloc(tsp->nnodes - 1, sizeof(double));
	char* row_name = (char*)calloc(100, sizeof(char));
	for (int h = 0; h < tsp->nnodes; h++) {
		int nnz = 0;
		for (int i = 0; i < tsp->nnodes; i++) {
			if (i == h)
				continue;
			vars[nnz] = xpos(i, h, tsp);
			coeffs[nnz] = 1;
			nnz++;
		}
		sprintf(row_name, "constraint(%d)", h + 1);
		int err;
		if ((err = CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &start_row_coefficients, vars, coeffs, NULL,
				      &row_name))) {
			printf("Error adding row: %d\n", err);
			return -1;
		}
	}
	free(row_name);
	free(coeffs);
	free(vars);

	// CPXwriteprob(env, lp, "model.lp", NULL);

	return 0;
}

int permutation_contains_value(struct tsp* tsp, int val, int limit)
{
	for (int i = 0; i < limit; i++) {
		if (tsp->solution_permutation[i] == val)
			return 1;
	}
	return 0;
}

int tsp_cplex_savesolution(struct tsp* tsp, CPXENVptr env, CPXLPptr lp)
{
	int res = 0;
	int ncols = CPXgetnumcols(env, lp);
	double* vars = malloc(sizeof(double) * ncols); // REVIEW this can get very large!

	if (CPXgetx(env, lp, vars, 0, ncols - 1)) {
		fprintf(stderr, "Error getting variable value\n");
		res = -1; // error
		goto out;
	}

	int current = 0;
	int i = 0;
	tsp->solution_permutation[current++] = i;
	while (current < tsp->nnodes) {
		for (int j = 0; j < tsp->nnodes; j++) {
			if (i != j) {
				int pos = xpos(i, j, tsp);
				int val = vars[pos] > 0.5 ? 1 : 0;
				if (val == 1 && !permutation_contains_value(tsp, j, current - 1) &&
				    !permutation_contains_value(tsp, i, current - 1)) {
					tsp->solution_permutation[current++] = j;
					i = j;
					break;
				}
			}
		}
	}

out:
	free(vars);
	return res;
}

#define EPS 1e-5
void tsp_cplex_buildsol(struct tsp* tsp, const double* xstar, int* succ, int* comp, int* ncomp)
{

#ifdef DEBUG
	int* degree = (int*)calloc(tsp->nnodes, sizeof(int));
	for (int i = 0; i < tsp->nnodes; i++) {
		for (int j = i + 1; j < tsp->nnodes; j++) {
			int k = xpos(i, j, tsp);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS)
				printf(" wrong xstar in build_sol()\n");
			if (xstar[k] > 0.5) {
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < tsp->nnodes; i++) {
		if (degree[i] != 2)
			printf("wrong degree in build_sol()\n");
	}
	free(degree);
#endif

	*ncomp = 0;
	for (int i = 0; i < tsp->nnodes; i++) {
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < tsp->nnodes; start++) {
		if (comp[start] >= 0)
			continue; // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done) // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < tsp->nnodes; j++) {
				if (i != j && xstar[xpos(i, j, tsp)] > 0.5 &&
				    comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start; // last arc to close the cycle

		// go to the next component
	}
}

double compute_rhs(int* comp, int len, int current_component)
{
	double rhs = -1;
	for (int i = 0; i < len; i++) {
		if (comp[i] == current_component) {
			rhs++;
		}
	}
	return rhs;
}

void tsp_cplex_addsec(struct tsp* tsp, CPXENVptr env, CPXLPptr lp, int ncomp, int* comp)
{
	int ncols = CPXgetnumcols(env, lp);

	int* index = malloc(sizeof(int) * ncols);
	double* value = malloc(sizeof(double) * ncols);
	char sense = 'L';
	char* cname = "";
	int izero = 0;
	for (int k = 1; k <= ncomp; k++) {
		int nnz = 0;
		double rhs = compute_rhs(comp, tsp->nnodes, k);
		for (int i = 0; i < tsp->nnodes; i++) {
			if (comp[i] != k)
				continue;
			for (int j = i + 1; j < tsp->nnodes; j++) {
				if (comp[j] != k)
					continue;
				index[nnz] = xpos(i, j, tsp);
				value[nnz] = 1.0;
				nnz++;
			}
		}
		int res = CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname);
		if (res) {
			printf("Can't add new row\n");
			goto free_buffers;
		}
	}

free_buffers:
	free(index);
	free(value);
}

void print_loops_file(struct tsp* tsp, int* succ, char* filename)
{
	int* visited = calloc(tsp->nnodes, sizeof(int));
	FILE* f = fopen(filename, "w");

	while (1) {
		int notvisited = -1;
		for (int i = 0; i < tsp->nnodes; i++) {
			if (!visited[i]) {
				notvisited = i;
				break;
			}
		}
		if (notvisited == -1)
			break;

		int current = notvisited;
		fprintf(f, "newloop\n");
		while (1) {
			visited[current] = 1;
			fprintf(f, "%lf, %lf\n", tsp->coords[current].x, tsp->coords[current].y);
			current = succ[current];
			if (current == notvisited) {
				fprintf(f, "%lf, %lf\n", tsp->coords[current].x, tsp->coords[current].y);
				break;
			}
		}
	}

	fclose(f);
	free(visited);
}

void tsp_cplex_patchonce(struct tsp* tsp, const int* succin, int* start, int ncomp, int* succout)
{
	// find the best patch
	double best_delta = 10e30;
	int besti = -1;
	int bestj = -1;
	int removed_comp = -1;
	for (int k1 = 1; k1 < ncomp; k1++) {
		for (int k2 = k1 + 1; k2 < ncomp + 1; k2++) {
			int current1 = start[k1];
			int current2 = start[k2];

			while (1) {
				while (1) {
					int i = current1;
					int j = current2;
					double c_isj = tsp->cost_matrix[flatten_coords(i, succin[j], tsp->nnodes)];
					double c_jsi = tsp->cost_matrix[flatten_coords(j, succin[i], tsp->nnodes)];
					double c_isi = tsp->cost_matrix[flatten_coords(i, succin[i], tsp->nnodes)];
					double c_jsj = tsp->cost_matrix[flatten_coords(j, succin[j], tsp->nnodes)];

					double delta = c_isj + c_jsi - c_isi - c_jsj;

					if (delta < best_delta) {
						best_delta = delta;
						besti = current1;
						bestj = current2;
						removed_comp = k2;
					}

					current2 = succin[current2];
					if (current2 == start[k2])
						break;
				}
				current1 = succin[current1];
				if (current1 == start[k1])
					break;
			}
		}
	}
	printf("patchonce: besti=%d, bestj=%d, removed=%d\n", besti, bestj, removed_comp);
	// execute the best patch
	memcpy(succout, succin, tsp->nnodes * sizeof(int));
	succout[succin[besti]] = succin[succin[bestj]];
	succout[succin[bestj]] = succin[succin[besti]];

	int temp = start[removed_comp];
	start[removed_comp] = start[ncomp];
	start[ncomp] = temp;
}

void tsp_compute_comp_start(struct tsp* tsp, int* comp, int ncomp, int* start)
{
	// initialize start
	for (int i = 0; i < ncomp + 1; i++)
		start[i] = -1;

	for (int i = 0; i < tsp->nnodes; i++) {
		int k = comp[i];
		if (start[k] == -1)
			start[k] = i;
	}
}

void tsp_cplex_patch_comp(struct tsp* tsp, const int* succin, int* comp, int ncomp, int* succout)
{
	int* start = malloc(sizeof(int) * (ncomp + 1));
	tsp_compute_comp_start(tsp, comp, ncomp, start);

	int* temp = malloc(sizeof(int) * tsp->nnodes);
	memcpy(temp, succin, tsp->nnodes * sizeof(int));
	printf("comp:\n");
	print_array_int(comp, tsp->nnodes);
	while (ncomp > 1) {
		printf("start:\n");
		print_array_int(start, ncomp + 1);
		printf("temp: \n");
		print_array_int(temp, tsp->nnodes);
		tsp_cplex_patchonce(tsp, temp, start, ncomp, succout);
		printf("succout: \n");
		print_array_int(succout, tsp->nnodes);
		memcpy(temp, succout, tsp->nnodes * sizeof(int));

		ncomp--;
		printf("ncomp = %d\n", ncomp);
	}

	free(temp);
	free(start);
}

int tsp_solve_cplex(struct tsp* tsp)
{
	tsp->solution_permutation = NULL;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	int res = 0;

	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) {
		printf("Error creating env: %d\n", error);
		res = -1;
		goto free_cplex;
	}
	CPXLPptr lp = CPXcreateprob(env, &error, "tsp");
	if (error) {
		printf("Error creating lp: %d\n", error);
		res = -1;
		goto free_prob;
	}

	if ((error = tsp_build_lpmodel(tsp, env, lp))) {
		printf("Erorr building model\n");
		res = -1;
		goto free_prob;
	}

#ifdef DEBUG
	CPXwriteprob(env, lp, "prob.lp", NULL);
#endif

	int* succ = malloc(sizeof(int) * tsp->nnodes);
	int* comp = malloc(sizeof(int) * tsp->nnodes);
	int ncomp;

	tsp_starttimer(tsp);

	int it = 0;

	while (1) {
		CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_getremainingseconds(tsp));
		it++;
		if (tsp_shouldstop(tsp)) {
			res = -1;
			goto timelimit_reached;
		}
#ifdef DEBUG
		int row = CPXgetnumrows(env, lp);
		fprintf(stderr, "Rows are %d\n", row);
#endif
		if ((error = CPXmipopt(env, lp))) {
			printf("Error while optimizing: %d\n", error);
			res = -1;
			goto free_buffers;
		}
		if (tsp_shouldstop(tsp)) {
			// this needs to be called again because the timelimit
			// could be reached during the execution of cplex.
			res = -1;
			goto timelimit_reached;
		}
		int ncols = CPXgetnumcols(env, lp);
		double* xvars = malloc(sizeof(double) * ncols);
		if (CPXgetx(env, lp, xvars, 0, ncols - 1)) {
			printf("Can't get vars\n");
			break;
		}
		tsp_cplex_buildsol(tsp, xvars, succ, comp, &ncomp);
		free(xvars);

		fprintf(stderr, "ncomp=%d\n", ncomp);

		if (ncomp == 1)
			break;
		tsp_cplex_addsec(tsp, env, lp, ncomp, comp);
#ifdef DEBUG
		char probname[50];
		sprintf(probname, DEBUGOUT_LPPROB, it);
		CPXwriteprob(env, lp, probname, NULL);
		sprintf(probname, DEBUGOUT_PARTIAL, it);
		// print_loops_file(tsp, succ, probname);
#endif

		// compute patching in case next
		// iteration doesn't have time to finish

		int* patched = malloc(sizeof(int) * tsp->nnodes);
		tsp_cplex_patch_comp(tsp, succ, comp, ncomp, patched);
		if (tsp_allocate_solution(tsp)) {
			fprintf(stderr, "Can't allocate solution\n");
		}
		if (tsp_succ_to_perm(tsp, patched, tsp->solution_permutation)) {
			fprintf(stderr, "Can't convert solution\n");
		}
#ifdef DEBUG
		sprintf(probname, DEBUGOUT_PATCHED, it);
		// print_loops_file(tsp, patched, probname);
#endif
		double cost = tsp_recompute_solution_arg(tsp, tsp->solution_permutation);
		fprintf(stderr, "Patched solution cost is %lf\n", cost);
		tsp->solution_value = cost;
		free(patched);
	}

	double objval;
	CPXgetobjval(env, lp, &objval);
	tsp->solution_value = objval;

	if (tsp_allocate_solution(tsp))
		res = -1;

	if (tsp_cplex_savesolution(tsp, env, lp)) {
		fprintf(stderr, "Can't get solution of lp\n");
		res = -1;
	}

timelimit_reached:
free_buffers:
	free(succ);
	free(comp);
free_prob:
	CPXfreeprob(env, &lp);
free_cplex:
	CPXcloseCPLEX(&env);
	return res;
}
