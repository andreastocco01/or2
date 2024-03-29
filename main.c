#include "eventlog.h"
#include "tsp.h"
#include "tsp_cplex.h"
#include "tsp_tabu.h"
#include "tsp_vns.h"
#include "util.h"
#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

void summary_and_exit(int signal);

time_t startTime;
int parseFriendly = 0;
int configPlot = 0;
int timeLimit = 0;
int plotCurrent = 0;
int runConfig = -1;
int childpid;
char* logfile = "log.txt";
struct tsp tsp;

int plot_instance(struct tsp* tsp)
{
	FILE* gpprocess = popen("gnuplot --persist", "w");
	if (gpprocess == NULL)
		return -1;

	fprintf(gpprocess, "$data << EOD\n");

	if (tsp->solution_permutation) {
		for (int i = 0; i < tsp->nnodes; i++) {
			fprintf(gpprocess, "%lf %lf\n", tsp->coords[tsp->solution_permutation[i]].x,
				tsp->coords[tsp->solution_permutation[i]].y);
		}
		fprintf(gpprocess, "%lf %lf\n", tsp->coords[tsp->solution_permutation[0]].x,
			tsp->coords[tsp->solution_permutation[0]].y);

		fprintf(gpprocess, "EOD\n");
		fprintf(gpprocess, "plot $data using 1:2 title \"dataset\" pt 7 ps 2 with "
				   "points, $data using 1:2 title \"solution\" with lines\n");
	} else {
		for (int i = 0; i < tsp->nnodes; i++) {
			fprintf(gpprocess, "%lf %lf\n", tsp->coords[i].x, tsp->coords[i].y);
		}
		fprintf(gpprocess, "EOD\n");
		fprintf(gpprocess, "plot $data using 1:2 pt 7 ps 2 with points\n");
	}

	fclose(gpprocess);
	return 0;
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
				tsp->edge_weight_type = (char*)malloc(sizeof(char) * strlen(value) + 1);
				strcpy(tsp->edge_weight_type, value);
			} else if (!strcmp(name, "TYPE") && strcmp(value, "TSP")) {
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
					fscanf(file, "%d %lf %lf\n", &index, &x, &y);
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
		} else if (!strcmp(argv[i], "--timelimit") || !strcmp(argv[i], "-t")) {
			timeLimit = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--plotcurrent") || !strcmp(argv[i], "-c")) {
			plotCurrent = 1;
		} else if (!strcmp(argv[i], "--parsefriendly")) {
			parseFriendly = 1;
		} else if (!strcmp(argv[i], "--config")) {
			runConfig = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--logfile")) {
			logfile = argv[++i];
		}
	}
}

void main_compute(int argc, char** argv)
{
	tsp_init(&tsp);
	if (eventlog_initialize(logfile)) {
		printf("Can't initialize logger\n");
		exit(-1);
	}

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

	/* tsp_compute_costs(&tsp, tsp_costfunction_euclidian); */
	tsp_compute_costs(&tsp, tsp_costfunction_att);

	if (runConfig == 0) {
		if (tsp_solve_vns(&tsp)) {
			perror("Can't solve greedy\n");
		}
	}
	if (runConfig == 1) {
		if (tsp_solve_tabu(&tsp, tenure_fixed)) {
			perror("Can't solve tabu\n");
		}
	}
	if (runConfig == 2) {
		if (tsp_solve_tabu(&tsp, tenure_sin)) {
			perror("Can't solve tabu\n");
		}
	}
	if (runConfig == 3) {
		if (tsp_solve_cplex(&tsp)) {
			perror("Can't solve cplex\n");
		}
	}

	fprintf(stderr, "Please specify a config to solve the problem\n");
	eventlog_close();
	exit(0);

	summary_and_exit(-1);
}

void parse_friendly_output()
{
	time_t currentTime = clock();
	double totalTime = (double)(currentTime - startTime) / CLOCKS_PER_SEC;

	printf("%lf;%lf\n", totalTime, tsp.solution_value);
}

void summary_and_exit(int signal)
{
	eventlog_close();

	double res;
	if (!tsp_check_solution(&tsp, &res)) {
		printf("The computed solution is invalid!\n");
	}

	if (parseFriendly) {
		parse_friendly_output();
		goto summary_finish;
	}

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
	}
summary_finish:
	exit(0);
}

void redirect_to_child(int signal)
{
	kill(childpid, signal);
}

int main(int argc, char** argv)
{
	parse_arguments(argc, argv);

	startTime = clock();
	childpid = fork();
	fprintf(stderr, "childpid = %d\n", childpid);
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
		if (!parseFriendly)
			printf("Child returned %s\n", strerror(res));
	}
	return 0;
}
