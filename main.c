#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "eventlog.h"
#include "tsp.h"
#include "tsp_instance.h"
#include "tsp_vns.h"

struct experiment_args {
	int parse_friendly;
	int do_plot;
	int runconfiguration;
	char* logfile;
};

struct experiment_args parse_arguments(int argc, char** argv)
{
	struct experiment_args args = {
	    .do_plot = 0,
	    .logfile = "log.txt",
	    .parse_friendly = 0,
	    .runconfiguration = -1
	};

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--plot") || !strcmp(argv[i], "-p")) {
			args.do_plot = 1;
		}  else if (!strcmp(argv[i], "--parsefriendly")) {
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

	fclose(gpprocess);
	return 0;
}

int run_experiment(struct tsp* tsp, int config)
{
    // add more configurations here
    return tsp_solve_vns(tsp);
	/* return tsp_solve_cplex(tsp); */
}

void print_parse_friendly_output(struct tsp* tsp)
{
	time_t currentTime = clock();
	double totalTime = (double)(currentTime - tsp->start_time) / CLOCKS_PER_SEC;

	printf("%lf;%lf\n", totalTime, tsp->solution_value);
}

int conclude_experiment(struct tsp* tsp, int is_parse_friendly, int do_plot)
{
    double res;
	if (!tsp_check_solution(tsp, &res)) {
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

int main(int argc, char** argv)
{
    struct experiment_args args = parse_arguments(argc, argv);
    struct tsp tsp;

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

    if(run_experiment(&tsp, args.runconfiguration)) {
        exit(-1);
    }

	conclude_experiment(&tsp, args.parse_friendly, args.do_plot);

    eventlog_close();

	return 0;
}
