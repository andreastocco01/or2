#include "tsp_tabu.h"

#include "tsp.h"
#include "tsp_greedy.h"
#include <bits/types/sigset_t.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int tenure_fixed(int nnodes, int iteration)
{
	int computed = nnodes / TENURE_FIXED_DIVISOR;
	return computed > TENURE_MIN ? computed : TENURE_MIN;
}

int tenure_sin(int nnodes, int iteration)
{
	int computed = sin(iteration * TENURE_SIN_SCALE) / TENURE_SIN_DIVISOR;
	return computed > TENURE_MIN ? computed : TENURE_MIN;
}

int tsp_solve_tabu(struct tsp* tsp, tsp_tenure tenure)
{
	if (tsp_allocate_solution(tsp))
		return -1;

	if (tsp_compute_costs(tsp))
		return -1;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	int* tabu_iteration = (int*)malloc(tsp->nnodes * sizeof(int));
	for (int i = 0; i < tsp->nnodes; i++)
		tabu_iteration[i] = -1;

	// starting solution is solved with a greedy approach
	int starting_node = rand() % tsp->nnodes;

	printf("starting from node %d\n", starting_node);
	tsp_solve_greedy(tsp, starting_node, tsp->solution_permutation, &tsp->solution_value);

	int* current_solution = malloc(tsp->nnodes * sizeof(int));
	double current_solution_value = tsp->solution_value;

	memcpy(current_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	tsp_add_current(tsp, current_solution_value);
	tsp_add_incumbent(tsp, current_solution_value);

	int current_iteration = 0;

	while (1) {
		for (int i = 0; i < tsp->nnodes - 2; i++) {
			double best_delta = -10e30;
			int best_i, best_j;
			for (int j = i + 2; j < tsp->nnodes; j++) {
				double delta = compute_delta(tsp, current_solution, i, j);
				if (delta > best_delta) {
					best_delta = delta;
					best_i = i;
					best_j = j;
				}
			}
			if (best_delta > 0) {
				// normal 2-opt
				tsp_2opt_swap(best_i + 1, best_j, current_solution);
				current_solution_value -= best_delta;
				tsp_add_current(tsp, current_solution_value);
				if (current_solution_value < tsp->solution_value) {
					tsp->solution_value = current_solution_value;
					tsp_add_incumbent(tsp, tsp->solution_value);

					memcpy(tsp->solution_permutation, current_solution, sizeof(int) * tsp->nnodes);
				}
			} else {
				if (tabu_iteration[i] != -1 &&
				    (current_iteration - tabu_iteration[i]) < tenure(tsp->nnodes, current_iteration)) {
					// the node is tabu. We need to skip
					continue;
				}
				tsp_2opt_swap(best_i + 1, best_j, current_solution);
				current_solution_value -= best_delta;
				// add to tabu list
				tabu_iteration[i] = current_iteration;
				tsp_add_current(tsp, current_solution_value);
			}
		}
		current_iteration++;
	}

	free(tabu_iteration);
	free(current_solution);
	return 0;
}
