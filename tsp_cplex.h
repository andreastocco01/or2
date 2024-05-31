#ifndef TSP_CPLEX_
#define TSP_CPLEX_

#include "tsp.h"
#include <ilcplex/cplex.h>

int tsp_build_lpmodel(struct tsp* tsp, CPXENVptr env, CPXLPptr lp);
int tsp_cplex_savesolution(struct tsp* tsp, CPXENVptr env, CPXLPptr lp);
int tsp_solve_benders(struct tsp* tsp, int patching);
int tsp_solve_branchcut(struct tsp* tsp, int warmstart, int fraccut, int post_heuristic);
int tsp_solve_branchcut_for_matheuristic(struct tsp* tsp,
					 CPXENVptr env,
					 CPXLPptr lp,
					 int warmstart,
					 int fraccut,
					 int post_heuristic);

/**
 * Converts a solution from the "permutation" format to the cplex format
 *
 * The array should be preallocated
 *
 * Returns 0 if success, 1 if fails.
 * */
int tsp_perm_to_cplex(const struct tsp* tsp, const int* perm, double* cplex_sol, int ncols);
int cplex_warm_start(struct tsp* tsp, CPXENVptr env, CPXLPptr lp);
int cplex_add_start(CPXENVptr env, CPXLPptr lp, double* solution, int ncols);

#endif
