#include "eventlog.h"
#include "tsp.h"
#include "tsp_cplex.h"
#include "tsp_diving.h"
#include "tsp_greedy.h"
#include "tsp_instance.h"
#include "tsp_localbranching.h"
#include "tsp_tabu.h"
#include "tsp_vns.h"
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct experiment_args {
	int parse_friendly;
	int do_plot;
	int runconfiguration;
	char* logfile;
};

struct experiment_args parse_arguments(int argc, char** argv)
{
	struct experiment_args args = {.do_plot = 0, .logfile = "log.txt", .parse_friendly = 0, .runconfiguration = -1};

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--plot") || !strcmp(argv[i], "-p")) {
			args.do_plot = 1;
		} else if (!strcmp(argv[i], "--parsefriendly")) {
			args.parse_friendly = 1;
		} else if (!strcmp(argv[i], "--config")) {
			args.runconfiguration = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--logfile")) {
			args.logfile = argv[++i];
		}
	}

	return args;
}

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

	pclose(gpprocess);
	return 0;
}

int run_experiment(struct tsp* tsp, int config)
{
	// multigreedy
	if (config == 0) {
		return tsp_solve_multigreedy(tsp, 0); // multigreedy
	}
	if (config == 1) {
		return tsp_solve_multigreedy(tsp, 1); // multigreedy + 2opt
	}

	// vns
	if (config == 200) {
		vns_setrange(1, 8);
		return tsp_solve_vns(tsp);
	}
	if (config == 201) {
		vns_setrange(1, 4);
		return tsp_solve_vns(tsp);
	}
	if (config == 202) {
		vns_setrange(1, 2);
		return tsp_solve_vns(tsp);
	}
	if (config == 203) {
		vns_setrange(2, 4);
		return tsp_solve_vns(tsp);
	}
	if (config == 204) {
		vns_setrange(4, 8);
		return tsp_solve_vns(tsp);
	}
	if (config == 205) {
		vns_setrange(4, 6);
		return tsp_solve_vns(tsp);
	}
	if (config == 206) {
		vns_setrange(6, 8);
		return tsp_solve_vns(tsp);
	}

	// tabu fixed tenure
	if (config == 3) {
		tenure_fixed_setdivisor(2.4);
		return tsp_solve_tabu(tsp, tenure_fixed);
	}
	if (config == 4) {
		tenure_fixed_setdivisor(2.2);
		return tsp_solve_tabu(tsp, tenure_fixed);
	}
	if (config == 5) {
		tenure_fixed_setdivisor(2.0);
		return tsp_solve_tabu(tsp, tenure_fixed);
	}
	if (config == 6) {
		tenure_fixed_setdivisor(1.8);
		return tsp_solve_tabu(tsp, tenure_fixed);
	}
	if (config == 7) {
		tenure_fixed_setdivisor(1.6);
		return tsp_solve_tabu(tsp, tenure_fixed);
	}

	// tabu sin tenure
	if (config == 8) {
		tenure_sin_setdivisor(200);
		tenure_sin_setscale(10);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 9) {
		tenure_sin_setdivisor(200);
		tenure_sin_setscale(20);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 10) {
		tenure_sin_setdivisor(100);
		tenure_sin_setscale(10);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 11) {
		tenure_sin_setdivisor(100);
		tenure_sin_setscale(20);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 12) {
		tenure_sin_setdivisor(50);
		tenure_sin_setscale(10);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 13) {
		tenure_sin_setdivisor(50);
		tenure_sin_setscale(20);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 14) {
		tenure_sin_setdivisor(25);
		tenure_sin_setscale(10);
		return tsp_solve_tabu(tsp, tenure_sin);
	}
	if (config == 15) {
		tenure_sin_setdivisor(25);
		tenure_sin_setscale(20);
		return tsp_solve_tabu(tsp, tenure_sin);
	}

	// benders
	if (config == 16) {
		return tsp_solve_benders(tsp, 1);
	}

	// branch and cut
	if (config == 17) {
		return tsp_solve_branchcut(tsp, 0, 0, 0);
	}
	if (config == 18) {
		return tsp_solve_branchcut(tsp, 1, 0, 0);
	}
	if (config == 19) {
		return tsp_solve_branchcut(tsp, 1, 1, 0);
	}
	if (config == 20) {
		return tsp_solve_branchcut(tsp, 0, 1, 0);
	}
	if (config == 21) {
		return tsp_solve_branchcut(tsp, 0, 0, 1);
	}
	if (config == 22) {
		return tsp_solve_branchcut(tsp, 1, 0, 1);
	}
	if (config == 23) {
		return tsp_solve_branchcut(tsp, 1, 1, 1);
	}
	if (config == 24) {
		return tsp_solve_branchcut(tsp, 0, 1, 1);
	}

	// diving
	if (config == 25) {
		return tsp_solve_diving(tsp, 0.9);
	}
	if (config == 26) {
		return tsp_solve_diving(tsp, 0.8);
	}
	if (config == 27) {
		return tsp_solve_diving(tsp, 0.7);
	}
	if (config == 28) {
		return tsp_solve_diving(tsp, 0.6);
	}
	if (config == 29) {
		return tsp_solve_diving(tsp, 0.5);
	}
	if (config == 30) {
		return tsp_solve_diving(tsp, 0.4);
	}

	// local branching
	if (config == 31) {
		return tsp_solve_localbranching(tsp, 10, 5);
	}
	if (config == 32) {
		return tsp_solve_localbranching(tsp, 20, 5);
	}
	return -1;
}

void print_parse_friendly_output(struct tsp* tsp)
{
	double elapsedseconds = tsp_getelapsedseconds(tsp);

	double sol;
	if (tsp->solution_permutation) {
		sol = tsp->solution_value;
	} else {
		sol = 10e30;
	}
	printf("%lf;%lf\n", elapsedseconds, sol);
}

int conclude_experiment(struct tsp* tsp, int is_parse_friendly, int do_plot)
{
	if (!tsp_is_solution(tsp)) {
		printf("The computed solution is not a solution!\n");
		return -1;
	}

	if (tsp->solution_permutation && !tsp_check_solution(tsp, NULL)) {
		printf("The computed solution is invalid!\n");
		return -1;
	}

	if (is_parse_friendly) {
		print_parse_friendly_output(tsp);
	} else {
		printf("Best solution: %lf\n", tsp->solution_value);
	}

	if (do_plot) {
		if (plot_instance(tsp)) {
			perror("Can't plot solution\n");
			return -1;
		}
	}

	return 0;
}

struct tsp* signal_tsp;

void handle_sigint(int signal)
{
	signal_tsp->force_stop = 1;
}

void setup_signals(struct tsp* tsp)
{
	signal_tsp = tsp;
	signal(SIGINT, handle_sigint);
}

int main(int argc, char** argv)
{
	struct experiment_args args = parse_arguments(argc, argv);
	struct tsp tsp;

	setup_signals(&tsp);

	tsp_init(&tsp);
	if (eventlog_initialize(args.logfile)) {
		printf("Can't initialize logger\n");
		exit(-1);
	}

	if (tsp_parse_arguments(argc, argv, &tsp)) {
		exit(-1);
	}

	if (tsp.model_source == 1) {
		tsp_loadinstance_random(&tsp);
	} else if (tsp.model_source == 2) {
		if (tsp_loadinstance_tsplib(&tsp) == -1)
			exit(-1);
	}

	/* tsp_compute_costs(&tsp, tsp_costfunction_euclidian); */
	tsp_compute_costs(&tsp, tsp_costfunction_att);

	if (run_experiment(&tsp, args.runconfiguration)) {
		fprintf(stderr, "Unable to find a solution\n");
	}

	conclude_experiment(&tsp, args.parse_friendly, args.do_plot);

	eventlog_close();

	tsp_free(&tsp);

	return 0;
}
