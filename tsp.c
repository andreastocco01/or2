#include "tsp.h"
#include <stdio.h>

void tsp_free(struct tsp* tsp)
{
	if (tsp->coords)
		free(tsp->coords);
}

void tsp_init(struct tsp* tsp)
{
	tsp->input_file = NULL;
	tsp->model_source = 0;
	tsp->coords = NULL;
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

void debug_print_coords(struct tsp* tsp)
{
	printf("%d nodes:\n", tsp->nnodes);
	for (int i = 0; i < tsp->nnodes; i++) {
		printf("%d: (%lf, %lf)\n", i + 1, tsp->coords[i].x,
		       tsp->coords[i].y);
	}
}
