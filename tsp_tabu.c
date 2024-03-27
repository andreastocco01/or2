#include "tsp_tabu.h"

#include "tsp.h"
#include "tsp_greedy.h"
#include <bits/types/sigset_t.h>
#include <math.h>
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
	int computed = sin(((double)iteration) / TENURE_SIN_DIVISOR) * TENURE_SIN_SCALE + TENURE_SIN_SCALE;
	return computed > TENURE_MIN ? computed : TENURE_MIN;
}

int is_tabu(int* tabu_iteration, int node, int current_iteration, int tenure)
{
	if (tabu_iteration[node] != -1 && ((current_iteration - tabu_iteration[node]) < tenure)) {
		return 1;
	}
	return 0;
}

double tsp_2opt_findbestswap_no_tabu(struct tsp* tsp,
				     int* solution,
				     int* best_i,
				     int* best_j,
				     int* tabu_iteration,
				     int tenure,
				     int current_iteration)
{
	double best_delta = -10e30;
	for (int i = 0; i < tsp->nnodes - 2; i++) {
		for (int j = i + 2; j < tsp->nnodes; j++) {
			double delta = compute_delta(tsp, solution, i, j);
			if (delta > best_delta) {
				// i have to choose the best nodes that aren't tabu!!!
				if (is_tabu(tabu_iteration, i, current_iteration, tenure) ||
				    is_tabu(tabu_iteration, j, current_iteration, tenure)) {
					continue;
				}
				best_delta = delta;
				*best_i = i;
				*best_j = j;
			}
		}
	}
	return best_delta;
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

#ifdef DEBUG
	fprintf(stderr, "starting from node %d\n", starting_node);
#endif
	tsp_solve_greedy(tsp, starting_node, tsp->solution_permutation, &tsp->solution_value);

	int* current_solution = malloc(tsp->nnodes * sizeof(int));
	double current_solution_value = tsp->solution_value;

	memcpy(current_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	tsp_add_current(tsp, current_solution_value);
	tsp_add_incumbent(tsp, current_solution_value);

	int current_iteration = 0;
	int ten;

	while (1) {
		// Intensification phase
		int best_i, best_j;
		double best_delta;
		while (1) {
			current_iteration++;
			ten = tenure(tsp->nnodes, current_iteration);
			best_delta = tsp_2opt_findbestswap_no_tabu(tsp, current_solution, &best_i, &best_j,
								   tabu_iteration, ten, current_iteration);
			if (best_delta <= 0)
				break; // local minimum

			tsp_2opt_solution(tsp, current_solution, &current_solution_value, best_i, best_j, best_delta);
		}

		// Diversification phase
		ten = tenure(tsp->nnodes, current_iteration);
		/* printf("tenure(%d) = %d\n", current_iteration, ten); */
		if (is_tabu(tabu_iteration, best_i, current_iteration, ten)) {
			// the node is tabu. We need to skip
			continue;
		}

		tsp_2opt_swap(best_i + 1, best_j, current_solution);
		current_solution_value -= best_delta;
		// add to tabu list
		tabu_iteration[best_i] = current_iteration;
		tsp_add_current(tsp, current_solution_value);
	}

	free(tabu_iteration);
	free(current_solution);
	return 0;
}
