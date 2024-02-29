#include "tsp.h"
#include <stdio.h>

void tsp_free(struct tsp* tsp)
{
	free(tsp->xcoord);
	free(tsp->ycoord);
}

void debug_print(struct tsp* tsp)
{
	printf("nnodes: %d\n", tsp->nnodes);

	printf("model_source: ");
	if (tsp->model_source == 1)
		printf("random\n");
	else if (tsp->model_source == 2)
		printf("input_file\n");
	else
		printf("NOT SET\n");

	printf("seed: %d\n", tsp->seed);
	printf("input_file: %s\n", tsp->input_file);
}
