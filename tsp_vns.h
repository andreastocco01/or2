#ifndef TSP_VNS_
#define TSP_VNS_

#include "tsp.h"

/**
 * Solve a tsp instance using the vns method
 * */
int tsp_solve_vns(struct tsp* tsp);

/**
 * Set the interval from which the number of kicks is generated
 * */
void vns_setrange(int min, int max);


#endif
