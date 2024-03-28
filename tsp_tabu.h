#ifndef TSP_TABU_
#define TSP_TABU_

#include "tsp.h"

// TENURES
typedef int (*tsp_tenure)(int nnodes, int iteration);

#define TENURE_MIN           10

#define TENURE_FIXED_DIVISOR 10
int tenure_fixed(int nnodes, int iteration);

#define TENURE_SIN_DIVISOR 1000
#define TENURE_SIN_SCALE   100
int tenure_sin(int nnodes, int iteration);

/**
 * Solve a tsp instance using the tabu method
 * This function runs forever.
 * */
int tsp_solve_tabu(struct tsp* tsp, tsp_tenure tenure);
int is_tabu(int* tabu_iteration, int node, int current_iteration, int tenure);
double tsp_2opt_findbestswap_no_tabu(struct tsp* tsp,
				     int* solution,
				     int* best_i,
				     int* best_j,
				     int* tabu_iteration,
				     int tenure,
				     int current_iteration);

#endif
