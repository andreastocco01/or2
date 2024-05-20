#include "tsp_cplex.h"
#include "ilcplex/cplex.h"
#include "mincut.h"
#include "tsp.h"
#include "tsp_greedy.h"
#include "util.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int xpos(int i, int j, const struct tsp* tsp)
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

int tsp_perm_to_cplex(const struct tsp* tsp, const int* perm, double* cplex_sol, int ncols)
{
	for (int i = 0; i < ncols; i++)
		cplex_sol[i] = 0.0;

	for (int i = 0; i + 1 < tsp->nnodes; i++) {
		int j = i + 1;
		int pos = xpos(perm[i], perm[j], tsp);
		cplex_sol[pos] = 1.0;
	}

	cplex_sol[xpos(perm[0], perm[tsp->nnodes - 1], tsp)] = 1.0;

	return 0;
}

#define EPS 1e-5
int tsp_cplex_buildsol(const struct tsp* tsp, const double* xstar, int* succ, int* comp, int* ncomp)
{
#ifdef DEBUG
	int* degree = (int*)calloc(tsp->nnodes, sizeof(int));
	for (int i = 0; i < tsp->nnodes; i++) {
		for (int j = i + 1; j < tsp->nnodes; j++) {
			int k = xpos(i, j, tsp);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) {
				printf(" wrong xstar in build_sol()\n");
				return -1;
			}
			if (xstar[k] > 0.5) {
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < tsp->nnodes; i++) {
		if (degree[i] != 2) {
			printf("wrong degree in build_sol()\n");
			return -1;
		}
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
	return 0;
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

void tsp_print_loops_file(struct tsp* tsp, int* succ, char* filename)
{
	int* visited = calloc(tsp->nnodes, sizeof(int));
	FILE* f = fopen(filename, "w");
	if (f == NULL) {
		fprintf(stderr, "Can't create file\n");
		return;
	}

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

void tsp_cplex_patchonce(const struct tsp* tsp, const int* succin, int* start, int ncomp, int* succout)
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
	// execute the best patch
	memcpy(succout, succin, tsp->nnodes * sizeof(int));
	succout[besti] = succin[bestj];
	succout[bestj] = succin[besti];

	int temp = start[removed_comp];
	start[removed_comp] = start[ncomp];
	start[ncomp] = temp;
}

void tsp_compute_comp_start(const struct tsp* tsp, int* comp, int ncomp, int* start)
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

void tsp_cplex_patch_comp(const struct tsp* tsp, const int* succin, int* comp, int ncomp, int* succout)
{
	int* start = malloc(sizeof(int) * (ncomp + 1));
	tsp_compute_comp_start(tsp, comp, ncomp, start);

	int* temp = malloc(sizeof(int) * tsp->nnodes);
	memcpy(temp, succin, tsp->nnodes * sizeof(int));
	while (ncomp > 1) {
		tsp_cplex_patchonce(tsp, temp, start, ncomp, succout);
		memcpy(temp, succout, tsp->nnodes * sizeof(int));

		ncomp--;
	}

	free(temp);
	free(start);
}

int tsp_solve_benders(struct tsp* tsp, int patching)
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
			free(xvars);
			break;
		}
		tsp_cplex_buildsol(tsp, xvars, succ, comp, &ncomp);
		free(xvars);

#ifdef DEBUG
		fprintf(stderr, "ncomp=%d\n", ncomp);
#endif

		if (ncomp == 1)
			break;
		tsp_cplex_addsec(tsp, env, lp, ncomp, comp);
#ifdef DEBUG
		char probname[50];
		sprintf(probname, DEBUGOUT_LPPROB, it);
		CPXwriteprob(env, lp, probname, NULL);
		sprintf(probname, DEBUGOUT_PARTIAL, it);
		tsp_print_loops_file(tsp, succ, probname);
#endif

		// compute patching in case next
		// iteration doesn't have time to finish

		if (patching) {
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
			tsp_print_loops_file(tsp, patched, probname);
#endif
			double cost = tsp_recompute_solution_arg(tsp, tsp->solution_permutation);
#ifdef DEBUG
			fprintf(stderr, "Patched solution cost is %lf\n", cost);
#endif
			tsp->solution_value = cost;
#ifdef DEBUG
			if (!tsp_check_solution(tsp, NULL)) {
				printf("Discrepancy in solution cost\n");
				exit(0);
			}
#endif
			tsp_2opt_swap_arg(tsp, tsp->solution_permutation, &tsp->solution_value);
			free(patched);
#ifdef DEBUG
			cost = tsp_recompute_solution_arg(tsp, tsp->solution_permutation);
			fprintf(stderr, "2-opted solution cost is %lf\n", cost);
			if (!tsp_check_solution(tsp, NULL)) {
				printf("Discrepancy in solution cost\n");
				exit(0);
			}
			sprintf(probname, DEBUGOUT_PATCHED2OPT, it);
			tsp_print_perm_file(tsp, tsp->solution_permutation, probname);
#endif
		}
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

struct callback_generate_sec_params {
	const struct tsp* tsp;
	int ncols;
	int post_heuristic;
};

static int callback_generate_sec_addsec(CPXCALLBACKCONTEXTptr context,
					int ncols,
					const struct tsp* tsp,
					int ncomp,
					int* comp)
{
	int res = 0;

	int* index = malloc(sizeof(int) * ncols);
	double* value = malloc(sizeof(double) * ncols);
	char sense = 'L';
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
		int res = CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value);
		if (res) {
			fprintf(stderr, "Can't reject solution and add new cut\n");
			res = -1;
			goto free_buffers;
		}
	}

free_buffers:
	free(index);
	free(value);
	return res;
}

static int CPXPUBLIC callback_generate_sec(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	int res = 0;
	const struct callback_generate_sec_params* params = (const struct callback_generate_sec_params*)userhandle;
	const struct tsp* tsp = params->tsp;
#ifdef DEBUG
	fprintf(stderr, "ncols is %d\n", params->ncols);
#endif

	// allocate buffers
	double* xstar = malloc(sizeof(double) * params->ncols);
	int* succ = malloc(sizeof(int) * tsp->nnodes);
	int* comp = malloc(sizeof(int) * tsp->nnodes);
	int ncomp;
	int* patched = malloc(sizeof(int) * tsp->nnodes);
	int* perm = malloc(sizeof(int) * tsp->nnodes);
	double* cplex_solution = malloc(sizeof(double) * params->ncols);
	double objval = CPX_INFBOUND;

	// retreive candidate solution
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, params->ncols - 1, &objval)) {
		fprintf(stderr, "callbackgetcandidatepoint failed!\n");
		res = -1; // triggers an error 1006 in cplex.
		goto free_buffers;
	}
#ifdef DEBUG
	fprintf(stderr, "Value of candidate is %lf\n", objval);
#endif

	// build a solution in the "successive" format
	if (tsp_cplex_buildsol(tsp, xstar, succ, comp, &ncomp)) {
		fprintf(stderr, "Can't build solution\n");
		res = -1;
		goto free_buffers;
	}
#ifdef DEBUG
	fprintf(stderr, "found %d components\n", ncomp);
#endif
	if (ncomp == 1) {
		// the solution is feasible
		res = 0;
		goto free_buffers;
	}

	// generate cuts
	if (callback_generate_sec_addsec(context, params->ncols, tsp, ncomp, comp)) {
		printf("Can't add SECs\n");
		res = -1;
		goto free_buffers;
	}

	if (params->post_heuristic) {

		// patch the solution
		tsp_cplex_patch_comp(tsp, succ, comp, ncomp, patched);

		// convert from successors to permutation
		if (tsp_succ_to_perm(tsp, patched, perm)) {
			fprintf(stderr, "Can't convert solution\n");
			res = -1;
			goto free_buffers;
		}

		// perform 2opt
		double cost = tsp_recompute_solution_arg(tsp, perm);
		tsp_2opt_swap_arg(tsp, perm, &cost);

		// convert from permutation to cplex format
		if (tsp_perm_to_cplex(tsp, perm, cplex_solution, params->ncols)) {
			fprintf(stderr, "Failed to convert from perm to cplex\n");
			res = -1;
			goto free_buffers;
		}

		int* ind = malloc(sizeof(int) * params->ncols);
		for (int i = 0; i < params->ncols; i++) {
			ind[i] = i;
		}

		int err = CPXcallbackpostheursoln(context, params->ncols, ind, cplex_solution, cost,
						  CPXCALLBACKSOLUTION_NOCHECK);
		if (err) {
			fprintf(stderr, "Failed to add heuristic solution: %d\n", err);
			res = -1;
			goto free_buffers;
		}

#ifdef DEBUG
		printf("Heuristic solution added\n");
#endif
		// free(val);
		free(ind);
	}

free_buffers:
	free(cplex_solution);
	free(perm);
	free(patched);
	free(xstar);
	free(succ);
	free(comp);
	return res;
}

struct violated_cut_callback_data {
	const struct callback_generate_sec_params* params;
	CPXCALLBACKCONTEXTptr context;
};

int add_cut_by_members(const struct tsp* tsp, CPXCALLBACKCONTEXTptr context, int* members, int membercount)
{
	int res = 0;
#ifdef DEBUG
	/* printf("Add cut for: "); */
	/* for (int i = 0; i < membercount; i++) { */
	/* 	printf("%d-", members[i]); */
	/* } */
	/* printf("\n"); */
#endif
	int max_edges = membercount * membercount; // TODO this is TOO generous
	char sense = 'L';
	int izero = 0;
	int purgeable = CPX_USECUT_PURGE;
	/* int purgeable = CPX_USECUT_FORCE; */
	double rhs = membercount - 1.0;
	int nnz = 0;
	int* index = malloc(sizeof(int) * max_edges);
	double* value = malloc(sizeof(double) * max_edges);
	for (int i = 0; i < membercount; i++) {
		for (int j = i + 1; j < membercount; j++) {
			int node_i = members[i];
			int node_j = members[j];
			index[nnz] = xpos(node_i, node_j, tsp);
			value[nnz] = 1.0;
			nnz++;
		}
	}
	int valid_only_local = 0;
	int addres = CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable,
					    &valid_only_local);
	if (addres)
		res = -1;
	free(index);
	free(value);
	return res;
}

int violated_cut_callback(double cut_value, int number_nodes, int* members, void* userhandle)
{
	const struct violated_cut_callback_data* data = (const struct violated_cut_callback_data*)userhandle;
	const struct callback_generate_sec_params* params = data->params;
	const struct tsp* tsp = params->tsp;
	return add_cut_by_members(tsp, data->context, members, number_nodes);
}

static int CPXPUBLIC callback_fraccut(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	CPXINT nodeid;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEUID, &nodeid);
	if (nodeid % 10)
		return 0; // so the callback is executed 1/10 of the times

	int res = 0;
	const struct callback_generate_sec_params* params = (const struct callback_generate_sec_params*)userhandle;
	const struct tsp* tsp = params->tsp;

	// allocate buffers
	double* xstar = malloc(sizeof(double) * params->ncols);
	int* succ = malloc(sizeof(int) * tsp->nnodes);
	int* comp = malloc(sizeof(int) * tsp->nnodes);
	double objval = CPX_INFBOUND;

	// retreive candidate solution
	if (CPXcallbackgetrelaxationpoint(context, xstar, 0, params->ncols - 1, &objval)) {
		fprintf(stderr, "callbackgetrelaxationpoint failed!\n");
		res = -1; // triggers an error 1006 in cplex.
		goto free_buffers;
	}
#ifdef DEBUG
	fprintf(stderr, "Value of candidate is %lf\n", objval);
#endif

	// input parameters
	int ncount = tsp->nnodes;
	int ecount = params->ncols;
	int* elist = malloc(sizeof(int) * ecount *
			    2); // I need 2 integers for each edge, representing the nodes it connects
	// TODO this can be done once for all the callbacks since
	// it will be always the same!
	for (int i = 0; i < tsp->nnodes; i++) {
		for (int j = i + 1; j < tsp->nnodes; j++) {
			int pos = xpos(i, j, tsp);
			elist[pos * 2] = i;
			elist[pos * 2 + 1] = j;
		}
	}
	double* x = xstar; // Value of each edge
	// output parameters
	int ncomp;
	int* compscount;
	int* comps;
	int connected_components = CCcut_connect_components(ncount, ecount, elist, x, &ncomp, &compscount, &comps);
	if (connected_components) {
		printf("Find connected components failed\n");
		return -1;
	}

	/* printf("Connected components are %d\n", ncomp); */
	/* for(int i=0 ; i < ncomp ; i++) */
	/* 	printf("%d - %d\n", i, compscount[i]); */

#ifdef DEBUG
	fprintf(stderr, "ncomp = %d\n", ncomp);
#endif
	if (ncomp == 1) {
		struct violated_cut_callback_data data = {.context = context, .params = params};
		CCcut_violated_cuts(ncount, ecount, elist, xstar, 1.9, violated_cut_callback, &data);
	} else {
		int counter = 0;
		for (int i = 0; i < ncomp; i++) {
			/* From the concorde documentation: */
			/* comps will return the nodes in the components (it will be an */
			/* ncount array, with the first compscount[0] elements making up */
			/* the first component, etc.) */
			add_cut_by_members(tsp, context, comps + counter, compscount[i]);
			counter += compscount[i];
		}
	}

	free(elist);
	/* free(x); */

free_buffers:
	free(xstar);
	free(succ);
	free(comp);

	return res;
}

static int CPXPUBLIC callback_dispatch(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
		return callback_generate_sec(context, contextid, userhandle);
	}
	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
		return callback_fraccut(context, contextid, userhandle);
	}

	return -1;
}

int cplex_add_start(CPXENVptr env, CPXLPptr lp, double* solution, int ncols)
{
	int beg = 0;
	int* varindices = malloc(sizeof(int) * ncols);
	for (int i = 0; i < ncols; i++) {
		varindices[i] = i;
	}
	char* name = "Heuristic start";
	/* int effort = CPX_MIPSTART_NOCHECK; */
	int effort = CPX_MIPSTART_AUTO;
	int res = CPXaddmipstarts(env, lp, 1, ncols, &beg, varindices, solution, &effort, &name);
	free(varindices);
	return res;
}

int cplex_warm_start(struct tsp* tsp, CPXENVptr env, CPXLPptr lp)
{
	double total = tsp->timelimit_secs;
	tsp->timelimit_secs = total / 10;
	// warm start: find a solution using an heuristic and pass it to CPLEX

	int greedyres = tsp_solve_multigreedy(tsp);
	if (greedyres) {
		fprintf(stderr, "Can't generate heuristic\n");
		return 1;
	} else {
		fprintf(stderr, "Generated heuristic with cost %lf\n", tsp->solution_value);
	}

	int ncols = CPXgetnumcols(env, lp);
	double* warm_solution = malloc(sizeof(double) * ncols);
	tsp_perm_to_cplex(tsp, tsp->solution_permutation, warm_solution, ncols);
	int addwarmres = cplex_add_start(env, lp, warm_solution, ncols);
	if (addwarmres) {
		fprintf(stderr, "Can't add mip start\n");
	}
	free(warm_solution);
	tsp->timelimit_secs = total - tsp->timelimit_secs;
	fprintf(stderr, "setting timelimit to %f\n", tsp->timelimit_secs);
	return 0;
}

int tsp_solve_branchcut(struct tsp* tsp, int warmstart, int fraccut, int post_heuristic)
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

	// set cplex parameters
#ifdef DEBUG
	CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
#endif
	int ncols = CPXgetnumcols(env, lp);

	struct callback_generate_sec_params params = {
	    .tsp = tsp,
	    .ncols = ncols,
	    .post_heuristic = post_heuristic,
	};

	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (fraccut)
		contextid |= CPX_CALLBACKCONTEXT_RELAXATION;
	if (CPXcallbacksetfunc(env, lp, contextid, callback_dispatch, &params)) {
		fprintf(stderr, "Can't set callback function!\n");
		res = -1;
		goto free_prob;
	}

	if (warmstart) {
		if (cplex_warm_start(tsp, env, lp)) {
			res = -1;
			goto free_prob;
		}
	}

	tsp_starttimer(tsp);
	CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_getremainingseconds(tsp));

	if (CPXmipopt(env, lp)) {
		fprintf(stderr, "Failed to solve the problem\n");
		res = -1;
		goto free_prob;
	}

	// TODO check the status

	double objval;
	CPXgetobjval(env, lp, &objval);
	tsp->solution_value = objval;

	if (tsp_allocate_solution(tsp))
		res = -1;

	if (tsp_cplex_savesolution(tsp, env, lp)) {
		fprintf(stderr, "Can't get solution of lp\n");
		res = -1;
	}

free_prob:
	CPXfreeprob(env, &lp);
free_cplex:
	CPXcloseCPLEX(&env);
	return res;
}

int tsp_solve_branchcut_matheuristic(struct tsp* tsp,
				     CPXENVptr env,
				     CPXLPptr lp,
				     int warmstart,
				     int fraccut,
				     int post_heuristic)
{
	tsp->solution_permutation = NULL;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	int res = 0;

#ifdef DEBUG
	CPXwriteprob(env, lp, "prob.lp", NULL);
#endif

	// set cplex parameters
#ifdef DEBUG
	CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
#endif
	int ncols = CPXgetnumcols(env, lp);

	struct callback_generate_sec_params params = {
	    .tsp = tsp,
	    .ncols = ncols,
	    .post_heuristic = post_heuristic,
	};

	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (fraccut)
		contextid |= CPX_CALLBACKCONTEXT_RELAXATION;
	if (CPXcallbacksetfunc(env, lp, contextid, callback_dispatch, &params)) {
		fprintf(stderr, "Can't set callback function!\n");
		res = -1;
		goto end;
	}

	if (warmstart) {
		if (cplex_warm_start(tsp, env, lp)) {
			res = -1;
			goto end;
		}
	}

	if (CPXmipopt(env, lp)) {
		fprintf(stderr, "Failed to solve the problem\n");
		res = -1;
		goto end;
	}

	// TODO check the status

	double objval;
	CPXgetobjval(env, lp, &objval);
	tsp->solution_value = objval;

	if (tsp_allocate_solution(tsp))
		res = -1;

	if (tsp_cplex_savesolution(tsp, env, lp)) {
		fprintf(stderr, "Can't get solution of lp\n");
		res = -1;
	}

end:
	return res;
}
