#include "tsp_cplex.h"
#include "tsp.h"

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

int tsp_cplex_getsolution(struct tsp* tsp, CPXENVptr env, CPXLPptr lp)
{
	int res = 0;
	int ncols = CPXgetnumcols(env, lp);
	double* vars = malloc(sizeof(double) * ncols); // REVIEW this can get very large!

	if(CPXgetx(env, lp, vars, 0, ncols-1)) {
		fprintf(stderr, "Error getting variable value\n");
		res = -1; // error
		goto out;
	}

	for (int i = 0; i < tsp->nnodes - 1; i++) {
		for (int j = i + 1; j < tsp->nnodes; j++) {
			int pos = xpos(i, j, tsp);
			int val = vars[pos] > 0.5 ? 1 : 0;
			printf("%d->%d: %d\n", i, j, val);
		}
	}

	// TODO actually get and store the solution

out:
	free(vars);
	return res;
}

int tsp_solve_cplex(struct tsp* tsp)
{
	if (tsp_allocate_solution(tsp))
		return -1;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	tsp_starttimer(tsp);

	int res = 0;

	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) {
		printf("Error creating env: %d\n", error);
		res = -1;
		goto free_cplex;
	}
	CPXLPptr lp = CPXcreateprob(env, &error, "tps");
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

	if ((error = CPXmipopt(env, lp))) {
		printf("Error while optimizing: %d\n", error);
		res = -1;
		goto free_prob;
	}

	double objval;
	CPXgetobjval(env, lp, &objval);
	printf("Result = %lf\n", objval);

	if(tsp_cplex_getsolution(tsp, env, lp)) {
		printf("Can't get solution of lp\n");
		res = -1;
		goto free_prob;
	}

free_prob:
	CPXfreeprob(env, &lp);
free_cplex:
	CPXcloseCPLEX(&env);
	return res;
}
