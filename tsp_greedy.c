#include "tsp_greedy.h"
#include "eventlog.h"
#include "tsp.h"
#include "tsp_tabu.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/**
 * Output buffers have to be preallocated
 * */
int tsp_solve_greedy(struct tsp* tsp, int starting_node, int* output_solution, double* output_value)
{
	if (starting_node < 0 || starting_node >= tsp->nnodes)
		return -1;

	if (output_solution == NULL || output_value == NULL)
		return -1;

	int* current_solution = malloc(sizeof(int) * tsp->nnodes);

	for (int i = 1; i < tsp->nnodes; i++)
		current_solution[i] = i;

	current_solution[0] = starting_node;
	current_solution[starting_node] = 0;

	double cumulative_dist = 0;

	for (int i = 0; i < tsp->nnodes - 1; i++) {
		double min_dist = 1e30;
		int min_index = -1;

		for (int j = i + 1; j < tsp->nnodes; j++) {
			double dist = tsp->cost_matrix[flatten_coords(current_solution[i], current_solution[j],
								      tsp->nnodes)];
			if (dist < min_dist) {
				min_dist = dist;
				min_index = j;
			}
		}

		int temp = current_solution[i + 1];
		current_solution[i + 1] = current_solution[min_index];
		current_solution[min_index] = temp;

		cumulative_dist += min_dist;
	}

	double backarc = tsp->cost_matrix[flatten_coords(current_solution[tsp->nnodes - 1], current_solution[0],
							 tsp->nnodes)];

	cumulative_dist += backarc;

	memcpy(output_solution, current_solution, sizeof(int) * tsp->nnodes);
	*output_value = cumulative_dist;

	free(current_solution);
	return 0;
}

int tsp_solve_multigreedy(struct tsp* tsp)
{
	if (tsp_allocate_solution(tsp))
		return -1;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	tsp->solution_value = 10e30;
	int* current_solution = malloc(sizeof(int) * tsp->nnodes);
	double current_solution_value = 10e30;

	int* to_extract = malloc(sizeof(int) * tsp->nnodes);
	for (int i = 0; i < tsp->nnodes; i++) {
		to_extract[i] = i;
	}
	int t = 0;
	int current_iteration = 0;

	tsp_starttimer(tsp);

	for (int i = 0; i < tsp->nnodes; i++) {
		// compute starting node
		int r = rand() % (tsp->nnodes - t);
		int starting_node = to_extract[r];
		to_extract[r] = to_extract[tsp->nnodes - t - 1];
		to_extract[tsp->nnodes - t - 1] = starting_node;
		t++;
#ifdef DEBUG
		fprintf(stderr, "Starting node %d/%d\n", starting_node + 1, tsp->nnodes);
#endif

		if (tsp_solve_greedy(tsp, starting_node, current_solution, &current_solution_value)) {
			fprintf(stderr, "Can't solve greedy!\n");
			return -1;
		}
		if (current_solution_value < tsp->solution_value) {
			tsp_save_solution(tsp, current_solution, current_solution_value);
			eventlog_logdouble("new_incumbent", current_iteration, current_solution_value);
		}
		eventlog_logdouble("new_current", current_iteration, current_solution_value);
		while (1) {
			if (tsp_shouldstop(tsp))
				goto free_solution_buffers;
			current_iteration++;
			int best_i, best_j;
			double best_delta = tsp_2opt_findbestswap(tsp, current_solution, &best_i, &best_j);
			if (best_delta <= 0)
				break;
			int isnewbest = tsp_2opt_solution(tsp, current_solution, &current_solution_value, best_i,
							  best_j, best_delta);
			if (isnewbest)
				eventlog_logdouble("new_incumbent", current_iteration, current_solution_value);
			eventlog_logdouble("new_current", current_iteration, current_solution_value);
		}
	}

free_solution_buffers:

	free(current_solution);
	free(to_extract);

	return 0;
}
