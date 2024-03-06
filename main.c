#include "tsp.h"
#include "util.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int configPlot = 0;

void plot_instance(struct tsp* tsp)
{
	FILE* gpprocess = popen("gnuplot --persist", "w");
	fprintf(gpprocess, "plot '-' using 1:2 title \"\" pt 7 ps 2 with "
			   "points \n");
	for (int i = 0; i < tsp->nnodes; i++) {
		fprintf(gpprocess, "%lf %lf\n", tsp->coords[i].x,
			tsp->coords[i].y);
	}
	fclose(gpprocess);
}

int load_instance_file(struct tsp* tsp)
{
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
				tsp_allocate_buffers(tsp);
			} else if (!strcmp(name, "EDGE_WEIGHT_TYPE")) {
				tsp->edge_weight_type = (char*)malloc(
				    sizeof(char) * strlen(value) + 1);
				strcpy(tsp->edge_weight_type, value);
			} else if (!strcmp(name, "TYPE") &&
				   strcmp(value, "TSP")) {
				perror("Wrong format\n");
				return -1;
			}
		} else {
			if (!strcmp(buffer, "NODE_COORD_SECTION")) {
				// start reading coordinates
				for (int i = 0; i < dim; i++) {
					int index;
					double x;
					double y;
					fscanf(file, "%d %lf %lf\n", &index, &x,
					       &y);
					assert(index == i + 1);
					tsp->coords[i].x = x;
					tsp->coords[i].y = y;
				}
			}
		}
	}

	free(buffer);
	return 0;
}

int load_instance_random(struct tsp* tsp)
{
	srand(tsp->seed);

	tsp_allocate_buffers(tsp);

	for (int i = 0; i < tsp->nnodes; i++) {
		tsp->coords[i].x = random01() * RANDOM_MAX_X;
		tsp->coords[i].y = random01() * RANDOM_MAX_Y;
	}

	return 0;
}

/*
 * Parses command line arguments for configuring the
 * execution
 * */
void parse_arguments(int argc, char** argv)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-plot")) {
			configPlot = 1;
		}
	}
}

int main(int argc, char** argv)
{
	struct tsp tsp;

	tsp_init(&tsp);

	if (tsp_parse_arguments(argc, argv, &tsp)) {
		return -1;
	}

	parse_arguments(argc, argv);

	if (tsp.model_source == 1) {
		load_instance_random(&tsp);
	} else if (tsp.model_source == 2) {
		if (load_instance_file(&tsp) == -1)
			return -1;
	}

#ifdef DEBUG
	debug_print_coords(&tsp);
#endif
	if (configPlot)
		plot_instance(&tsp);

#ifdef DEBUG
	debug_print(&tsp);
#endif

	tsp_free(&tsp);
	return 0;
}
