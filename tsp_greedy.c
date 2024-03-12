#include "tsp_greedy.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Output buffers have to be preallocated
 * */
int tsp_solve_greedy(struct tsp* tsp,
		     int starting_node,
		     int* output_solution,
		     double* output_value)
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
			double dist = tsp->cost_matrix[flatten_coords(
			    current_solution[i], current_solution[j],
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

	double backarc = tsp->cost_matrix[flatten_coords(
	    current_solution[tsp->nnodes - 1], current_solution[0],
	    tsp->nnodes)];

	cumulative_dist += backarc;

	memcpy(output_solution, current_solution, sizeof(int) * tsp->nnodes);
	*output_value = cumulative_dist;

	free(current_solution);
	return 0;
}

int tsp_solve_multigreedy(struct tsp* tsp,
			  int* output_solution,
			  double* output_value)
{
	if (tsp->nnodes <= 0)
		return -1;

	int* best = malloc(sizeof(int) * tsp->nnodes);
	int* current = malloc(sizeof(int) * tsp->nnodes);
	double current_dist = 10e30;
	double best_dist = 10e30;

	for (int i = 0; i < tsp->nnodes; i++) {
#ifdef DEBUG
		printf("Starting node %d/%d\n", i + 1, tsp->nnodes);
#endif
		tsp_solve_greedy(tsp, i, current, &current_dist);
		tsp_2opt_solution(tsp, current, &current_dist);

		if (current_dist < best_dist) {
			best_dist = current_dist;
			memcpy(best, current, sizeof(int) * tsp->nnodes);
		}
	}

	memcpy(output_solution, best, sizeof(int) * tsp->nnodes);
	*output_value = best_dist;

	free(best);
	free(current);

	return 0;
}

int tsp_solve_multigreedy_save(struct tsp* tsp)
{
	if (tsp_allocate_solution(tsp))
		return -1;
	if (tsp_allocate_costs(tsp))
		return -1;
	if (tsp_compute_costs(tsp))
		return -1;

	if (tsp_solve_multigreedy(tsp, tsp->solution_permutation,
				  &tsp->solution_value))
		return -1;

	return 0;
}

int tsp_solve_greedy_save(struct tsp* tsp, int starting_node)
{
	if (tsp_allocate_solution(tsp))
		return -1;
	if (tsp_allocate_costs(tsp))
		return -1;
	if (tsp_compute_costs(tsp))
		return -1;

	if (tsp_solve_greedy(tsp, starting_node, tsp->solution_permutation,
			     &tsp->solution_value))
		return -1;

	return 0;
}
