#ifndef TSP_CPLEX_
#define TSP_CPLEX_

#include "tsp.h"
#include <ilcplex/cplex.h>

int tsp_build_lpmodel(struct tsp* tsp, CPXENVptr env, CPXLPptr lp);
int tsp_cplex_savesolution(struct tsp* tsp, CPXENVptr env, CPXLPptr lp);
int tsp_solve_benders(struct tsp* tsp, int patching);
int tsp_solve_branchcut(struct tsp* tsp);

#endif
