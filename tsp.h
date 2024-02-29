#ifndef TSP_H_
#define TSP_H_

#include <stdlib.h>

struct tsp {
	// instance data
	int nnodes;
	double* xcoord;
	double* ycoord;

	int model_source; // 1 if random, 2 if input_file

	// random generator data
	int seed;

	// input file data
	char* input_file;
};

void tsp_free(struct tsp* tsp);
void debug_print(struct tsp* tsp);

#endif // TSP_H_
