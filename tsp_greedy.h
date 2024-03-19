#ifndef TSP_GREEDY_
#define TSP_GREEDY_

#include "tsp.h"
#include <time.h>

int tsp_solve_greedy(struct tsp* tsp, int starting_node, int* output_solution, double* output_value);

int tsp_solve_multigreedy(struct tsp* tsp, int* output_solution, double* output_value);

int tsp_solve_multigreedy_init(struct tsp* tsp);

void tsp_2opt_swap(int left, int right, int* solution);

#endif
