#include "tsp_instance.h"
#include "util.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int tsp_loadinstance_random(struct tsp* tsp) {

	srand(tsp->seed);

	if (tsp_allocate_buffers(tsp) != 0)
		return -1;

	for (int i = 0; i < tsp->nnodes; i++) {
		tsp->coords[i].x = random01() * RANDOM_MAX_X;
		tsp->coords[i].y = random01() * RANDOM_MAX_Y;
	}

	return 0;
}

int tsp_loadinstance_tsplib(struct tsp* tsp)
{
	int res = 0;
	FILE* file = fopen(tsp->input_file, "r");
	if (file == NULL)
		return -1;
	char* buffer = malloc(64);

	size_t n = 64;
	size_t read;

	int dim = 0;

	while ((read = getline(&buffer, &n, file)) != -1) {
		buffer[read - 1] = 0; // remove the newline
		if (strstr(buffer, ":") != NULL) {
			// this is a property
			char* name = strtok(buffer, " : ");
			char* value = strtok(NULL, " : ");

			if (!strcmp(name, "DIMENSION")) {
				dim = atoi(value);
				tsp->nnodes = dim;
				if (tsp_allocate_buffers(tsp) != 0) {
					res = -1;
					goto free_buffer;
				}
			} else if (!strcmp(name, "EDGE_WEIGHT_TYPE")) {
				tsp->edge_weight_type = (char*)malloc(sizeof(char) * strlen(value) + 1);
				strcpy(tsp->edge_weight_type, value);
			} else if (!strcmp(name, "TYPE") && strcmp(value, "TSP")) {
				perror("Wrong format\n");
				res = -1;
				goto free_buffer;
			}
		} else {
			if (!strcmp(buffer, "NODE_COORD_SECTION")) {
				// start reading coordinates
				for (int i = 0; i < dim; i++) {
					int index;
					double x;
					double y;
					fscanf(file, "%d %lf %lf\n", &index, &x, &y);
					assert(index == i + 1);
					tsp->coords[i].x = x;
					tsp->coords[i].y = y;
				}
			}
		}
	}

free_buffer:
	free(buffer);
	return res;
}
