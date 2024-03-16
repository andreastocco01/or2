#ifndef TSP_GREEDY_C
#define TSP_GREEDY_C

#include "tsp.h"
#include <time.h>

int tsp_solve_greedy(struct tsp* tsp,
		     int starting_node,
		     int* output_solution,
		     double* output_value);

int tsp_solve_multigreedy(struct tsp* tsp,
			  int* output_solution,
			  double* output_value);

int tsp_solve_multigreedy_init(struct tsp* tsp);

int tsp_solve_greedy_save(struct tsp* tsp, int starting_node);

int tsp_solve_tabu(struct tsp* tsp);

#endif
