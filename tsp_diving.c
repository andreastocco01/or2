#include "tsp_diving.h"
#include "ilcplex/cplex.h"
#include "tsp.h"
#include "tsp_cplex.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>

int fix_edges(struct tsp* tsp, CPXENVptr env, CPXLPptr lp, double percentage, int* solution_permutation)
{
	int ncols = CPXgetnumcols(env, lp);
	double* cplex_solution = malloc(sizeof(double) * ncols);
	if (tsp_perm_to_cplex(tsp, solution_permutation, cplex_solution, ncols)) {
		fprintf(stderr, "Failed to convert from perm to cplex\n");
		free(cplex_solution);
		return 1;
	}
	srand(tsp->seed);
	char bound;
	double bound_value;
	for (int i = 0; i < ncols; i++) {
		double r = ((double)rand()) / RAND_MAX;
		if (r <= percentage) {
			if (cplex_solution[i] > 0.1) {
				bound = 'L';
				bound_value = 1.0;
				if (CPXchgbds(env, lp, 1, &i, &bound, &bound_value) != 0) {
					printf("Can't fix the bound\n");
					return 1;
				}
			} else {
				bound = 'U';
				bound_value = 0.0;
				if (CPXchgbds(env, lp, 1, &i, &bound, &bound_value) != 0) {
					printf("Can't fix the bound\n");
					return 1;
				}
			}
		}
	}
	free(cplex_solution);
	return 0;
}

int unfix_edges(struct tsp* tsp, CPXENVptr env, CPXLPptr lp)
{
	int ncols = CPXgetnumcols(env, lp);

	char bound;
	double bound_value;
	for (int i = 0; i < ncols; i++) {
		bound = 'L';
		bound_value = 0.0;
		if (CPXchgbds(env, lp, 1, &i, &bound, &bound_value) != 0) {
			printf("Can't fix the bound\n");
			return 1;
		}
		bound = 'U';
		bound_value = 1.0;
		if (CPXchgbds(env, lp, 1, &i, &bound, &bound_value) != 0) {
			printf("Can't fix the bound\n");
			return 1;
		}
	}
	return 0;
}

int tsp_solve_diving(struct tsp* tsp, double percentage)
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
	int* best_solution = malloc(sizeof(int) * tsp->nnodes);
	int* starting_solution = malloc(sizeof(int) * tsp->nnodes);
	memcpy(best_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	memcpy(starting_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	double best_obj = tsp->solution_value;

	while (1) {
		CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_getremainingseconds(tsp));
		if (tsp_shouldstop(tsp)) {
			res = -1;
			goto time_limit_reached;
		}

		// fix edges
		if (fix_edges(tsp, env, lp, percentage, starting_solution)) {
			printf("Unable to fix edges\n");
			res = -1;
			goto free_prob;
		}

		if (tsp_solve_branchcut_matheuristic(tsp, env, lp, 0, 1, 1) != 0) {
			printf("Unable to solve branch and cut\n");
			res = -1;
			goto free_prob;
		}

		if (tsp->solution_value < best_obj) {
			memcpy(best_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
			best_obj = tsp->solution_value;
		}

		// unfix edges
		if (unfix_edges(tsp, env, lp)) {
			printf("Unable to unfix edges\n");
			res = -1;
			goto free_prob;
		}
	}

time_limit_reached:
	memcpy(tsp->solution_permutation, best_solution, sizeof(int) * tsp->nnodes);
	tsp->solution_value = best_obj;
free_prob:
	CPXfreeprob(env, &lp);
free_cplex:
	CPXcloseCPLEX(&env);
	return res;
}
