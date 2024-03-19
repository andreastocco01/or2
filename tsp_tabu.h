#ifndef TSP_TABU_
#define TSP_TABU_

#include "tsp.h"

// TENURES
typedef int (*tsp_tenure)(int nnodes, int iteration);

#define TENURE_MIN           10

#define TENURE_FIXED_DIVISOR 10
int tenure_fixed(int nnodes, int iteration);

#define TENURE_SIN_DIVISOR 10
#define TENURE_SIN_SCALE   50
int tenure_sin(int nnodes, int iteration);

int tsp_solve_tabu(struct tsp* tsp, tsp_tenure tenure);

#endif
