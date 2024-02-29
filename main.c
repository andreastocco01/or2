#include "tsp.h"
#include "util.h"
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
	return 0;
}

int load_instance_random(struct tsp* tsp)
{
	srand(tsp->seed);

	tsp->coords = (struct point*)malloc(sizeof(struct point) * tsp->nnodes);

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

	debug_print(&tsp);

	if (tsp.model_source == 1) {
		load_instance_random(&tsp);
	} else if (tsp.model_source == 2) {
		load_instance_file(&tsp);
	}

	debug_print_coords(&tsp);
	if (configPlot)
		plot_instance(&tsp);

	tsp_free(&tsp);
	return 0;
}
