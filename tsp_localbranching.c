#include "tsp_localbranching.h"
#include "ilcplex/cplex.h"
#include "tsp.h"
#include "tsp_cplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int add_localbranching_constraint(const struct tsp* tsp, CPXENVptr env, CPXLPptr lp, int* solution_permutation, int k)
{
	int ncols = CPXgetnumcols(env, lp);
	double* cplex_solution = malloc(sizeof(double) * ncols);
	int* indexes = malloc(sizeof(int) * ncols);
	if (tsp_perm_to_cplex(tsp, solution_permutation, cplex_solution, ncols)) {
		fprintf(stderr, "Failed to convert from perm to cplex\n");
		free(cplex_solution);
		return 1;
	}

	for (int i = 0; i < ncols; i++) {
		if (cplex_solution[i] > 0.5)
			cplex_solution[i] = 1.0;
		else
			cplex_solution[i] = 0.0;

		indexes[i] = i;
	}

	char sense = 'G';

	double rhs = tsp->nnodes - k;

	char* constraint_name = "lbc";

	int zero = 0;

	int res = CPXaddrows(env, lp, 0, 1, ncols, &rhs, &sense, &zero, indexes, cplex_solution, NULL,
			     &constraint_name);

	free(indexes);
	free(cplex_solution);
	return res;
}

int remove_lastconstraint(CPXENVptr env, CPXLPptr lp)
{
	int nrows = CPXgetnumrows(env, lp);

	return CPXdelrows(env, lp, nrows - 1, nrows - 1);
}

int tsp_solve_localbranching(struct tsp* tsp, int k_start, int k_delta)
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

	// adding warm start to cplex
	if (cplex_warm_start(tsp, env, lp)) {
		res = -1;
		goto free_prob;
	}

	tsp_starttimer(tsp);

	int ncols = CPXgetnumcols(env, lp);
	int* best_solution = malloc(sizeof(int) * tsp->nnodes);
	double* cplex_best_solution = malloc(sizeof(double) * ncols);
	memcpy(best_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	double best_obj = tsp->solution_value;

	int current_k = k_start;

	while (1) {
		CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_getremainingseconds(tsp));
		if (tsp_shouldstop(tsp)) {
			goto time_limit_reached;
		}

		if (add_localbranching_constraint(tsp, env, lp, best_solution, current_k)) {
			fprintf(stderr, "Unable to add local branching constraint\n");
			res = -1;
			goto free_prob;
		}

		fprintf(stderr, "Solving b&c with k=%d\n", current_k);

		if (tsp_solve_branchcut_for_matheuristic(tsp, env, lp, 0, 1, 1) != 0) {
			fprintf(stderr, "Unable to solve branch and cut\n");
			res = -1;
			goto free_prob;
		}

		if (remove_lastconstraint(env, lp)) {
			fprintf(stderr, "Unable to remove last constraint\n");
			res = -1;
			goto free_prob;
		}

		int inc_improved = 0;
		if (tsp->solution_value < best_obj) {
			fprintf(stderr, "Found new incumbent = %lf\n", tsp->solution_value);
			memcpy(best_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
			best_obj = tsp->solution_value;
			tsp_perm_to_cplex(tsp, best_solution, cplex_best_solution, ncols);
			cplex_add_start(env, lp, cplex_best_solution, ncols);
			inc_improved = 1;
		}

		if (current_k > tsp->nnodes)
			goto optimal_found;

		if (!inc_improved)
			current_k += k_delta;
	}

optimal_found:
time_limit_reached:
	memcpy(tsp->solution_permutation, best_solution, sizeof(int) * tsp->nnodes);
	tsp->solution_value = best_obj;
	free(cplex_best_solution);
	free(best_solution);
free_prob:
	CPXfreeprob(env, &lp);
free_cplex:
	CPXcloseCPLEX(&env);
	return res;
}
