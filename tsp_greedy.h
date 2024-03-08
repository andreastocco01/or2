#ifndef TSP_GREEDY_C
#define TSP_GREEDY_C

#include "tsp.h"

int tsp_solve_greedy(struct tsp* tsp,
		     int starting_node,
		     int* output_solution,
		     double* output_value);


int tsp_solve_multigreedy(struct tsp* tsp,
			  int* output_solution,
			  double* output_value);


int tsp_solve_multigreedy_save(struct tsp* tsp);


int tsp_solve_greedy_save(struct tsp* tsp, int starting_node);

#endif
