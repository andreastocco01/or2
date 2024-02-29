#ifndef TSP_H_
#define TSP_H_

#include <stdlib.h>

struct tsp {
	int nnodes;
	double* xcoord;
	double* ycoord;
};

void tsp_free(struct tsp* tsp);

#endif // TSP_H_
