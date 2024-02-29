#ifndef TSP_H_
#define TSP_H_

#include <stdlib.h>

#define RANDOM_MAX_X 10000
#define RANDOM_MAX_Y 10000

struct point {
	double x;
	double y;
};

struct tsp {
	// instance data
	int nnodes;
	struct point* coords;

	int model_source; // 1 if random, 2 if input_file

	// random generator data
	int seed;

	// input file data
	char* input_file;
};

void tsp_init(struct tsp* tsp);
void tsp_free(struct tsp* tsp);
void debug_print(struct tsp* tsp);
void debug_print_coords(struct tsp* tsp);

#endif // TSP_H_
