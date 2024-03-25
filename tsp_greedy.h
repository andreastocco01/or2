#ifndef TSP_GREEDY_
#define TSP_GREEDY_

#include "tsp.h"
#include <time.h>

/**
 * Solve a tsp instance a pure greedy approach
 * */
int tsp_solve_greedy(struct tsp* tsp, int starting_node, int* output_solution, double* output_value);

/**
 * Solve a tsp instance using the greedy multistart approach and the 2opt
 * optimization
 *
 * TODO should this be moved to another file?
 * */
int tsp_solve_multigreedy(struct tsp* tsp, int* output_solution, double* output_value);

/**
 * Allocate buffers and solve with multigreedy
 * */
int tsp_solve_multigreedy_init(struct tsp* tsp);

#endif
