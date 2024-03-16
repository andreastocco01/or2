#include "tsp.h"
#include "tsp_greedy.h"
#include "util.h"
#include <assert.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

void summary_and_exit(int signal);

int configPlot = 0;
int timeLimit = 0;
int plot_current = 0;
int childpid;
struct tsp tsp;

int plot_instance(struct tsp* tsp)
{
	FILE* gpprocess = popen("gnuplot --persist", "w");
	if (gpprocess == NULL)
		return -1;

	fprintf(gpprocess, "$data << EOD\n");

	if (tsp->solution_permutation) {
		for (int i = 0; i < tsp->nnodes; i++) {
			fprintf(gpprocess, "%lf %lf\n",
				tsp->coords[tsp->solution_permutation[i]].x,
				tsp->coords[tsp->solution_permutation[i]].y);
		}
		fprintf(gpprocess, "%lf %lf\n",
			tsp->coords[tsp->solution_permutation[0]].x,
			tsp->coords[tsp->solution_permutation[0]].y);

		fprintf(gpprocess, "EOD\n");
		fprintf(
		    gpprocess,
		    "plot $data using 1:2 title \"dataset\" pt 7 ps 2 with "
		    "points, $data using 1:2 title \"solution\" with lines\n");
	} else {
		for (int i = 0; i < tsp->nnodes; i++) {
			fprintf(gpprocess, "%lf %lf\n", tsp->coords[i].x,
				tsp->coords[i].y);
		}
		fprintf(gpprocess, "EOD\n");
		fprintf(gpprocess,
			"plot $data using 1:2 pt 7 ps 2 with points\n");
	}

	fclose(gpprocess);
	return 0;
}

int plot_array(double* array, int size, const char* title, FILE* file)
{
	if (file == NULL)
		return -1;

	fprintf(file, "$data << EOD\n");

	for (int i = 0; i < size; i++) {
		fprintf(file, "%d %lf\n", i, array[i]);
	}

	fprintf(file, "EOD\n");
	fprintf(file,
		"plot $data using 1:2 title \"%s\" pt 7 ps 2 with lines\n",
		title);

	fclose(file);
	return 0;
}

int plot_incumbents(struct tsp* tsp)
{
	FILE* gpprocess = popen("gnuplot --persist", "w");
	return plot_array(tsp->incumbents, tsp->incumbent_next_index,
			  "incumbent", gpprocess);
}

int plot_current_solutions(struct tsp* tsp)
{
	FILE* gpprocess = popen("gnuplot --persist", "w");
	return plot_array(tsp->current_solutions,
			  tsp->current_solution_next_index, "current solution",
			  gpprocess);
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
				if (tsp_allocate_buffers(tsp) != 0) {
					return -1;
				}
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

	if (tsp_allocate_buffers(tsp) != 0)
		return -1;

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
		if (!strcmp(argv[i], "--plot") || !strcmp(argv[i], "-p")) {
			configPlot = 1;
		} else if (!strcmp(argv[i], "--timelimit") ||
			   !strcmp(argv[i], "-t")) {
			timeLimit = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--plotcurrent") ||
			   !strcmp(argv[i], "-c")) {
			plot_current = 1;
		}
	}
}

void main_compute(int argc, char** argv)
{
	tsp_init(&tsp);

	if (tsp_parse_arguments(argc, argv, &tsp)) {
		exit(-1);
	}

	if (tsp.model_source == 1) {
		load_instance_random(&tsp);
	} else if (tsp.model_source == 2) {
		if (load_instance_file(&tsp) == -1)
			exit(-1);
	}

#ifdef DEBUG
	debug_print(&tsp);
	debug_print_coords(&tsp);
#endif

	/* printf("------------GREEDY-----------\n"); */
	/* if (tsp_solve_multigreedy_init(&tsp)) { */
	/* 	perror("Can't solve greedy\n"); */
	/* } */
	printf("------------TABU-----------\n");
	if (tsp_solve_tabu(&tsp)) {
		perror("Can't solve tabu search\n");
	}
	summary_and_exit(-1);
}

void summary_and_exit(int signal)
{
	if (signal == -1) {
		kill(getppid(), SIGINT);
		printf("Execution terminated\n");
	} else if (signal == SIGUSR1) {
		printf("Timelimit reached\n");
	} else if (signal == SIGKILL) {
		printf("Execution terminated by the user\n");
	}

#ifdef DEBUG
	debug_print(&tsp);
#endif

	if (configPlot) {
		if (plot_instance(&tsp)) {
			perror("Can't plot solution\n");
		}
		if (plot_incumbents(&tsp)) {
			perror("Can't plot incumbents\n");
		}
		if (plot_current == 1 && plot_current_solutions(&tsp)) {
			perror("Can't plot incumbents\n");
		}
	}
	exit(0);
}

void redirect_to_child(int signal)
{
	kill(childpid, signal);
}

int main(int argc, char** argv)
{
	parse_arguments(argc, argv);

	childpid = fork();
	printf("childpid = %d\n", childpid);
	if (!childpid) {
		// ===== CHILD
		// signal for terminating with timelimit
		signal(SIGUSR1, summary_and_exit);
		// signal for terminating with ctrl+c
		signal(SIGINT, summary_and_exit);
		main_compute(argc, argv);
	} else {
		// ===== PARENT
		signal(SIGINT, redirect_to_child);
		if (timeLimit != 0) {
			sleep(timeLimit);
			int res = kill(childpid, SIGUSR1);
		}
		int res;
		wait(&res);
		printf("Child returned %s\n", strerror(res));
	}
	return 0;
}
