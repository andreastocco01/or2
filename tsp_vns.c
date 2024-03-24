#include "tsp_vns.h"
#include "tsp_greedy.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int compar(const void* a, const void* b)
{
	return *((int*)a) - *((int*)b);
}

void generate_3opt_positions(struct tsp* tsp, int* positions, int size)
{
	int count = 0;
	while (count < size) {
		int found = 0;
		int r = rand() % tsp->nnodes;

		for (int i = 0; i < count; i++) {
			if (positions[i] == r) {
				found = 1;
				break;
			}
		}

		if (found != 1) {
			positions[count++] = r;
		}
	}

	qsort(positions, size, sizeof(int), compar);
}

void tsp_3opt_swap(int i, int j, int k, int* current_solution, int* new_solution, int size)
{
	memcpy(new_solution, current_solution, sizeof(int) * i + 1); // elements from 0 to i

	int pos = i + 1;
	int z = j;
	// elements from j to i + 1
	while (z >= i + 1) {
		new_solution[pos++] = current_solution[z--];
	}

	z = k;
	while (z >= j + 1) {
		new_solution[pos++] = current_solution[z--];
	}

	memcpy(new_solution + pos, current_solution + k + 1, sizeof(int) * (size - k - 1));
}

int tsp_solve_vns(struct tsp* tsp)
{
	if (tsp_allocate_solution(tsp))
		return -1;

	if (tsp_compute_costs(tsp))
		return -1;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

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
				int positions[3];
				generate_3opt_positions(tsp, positions, 3); // i -> j -> k
			}
		}
		current_iteration++;
	}

	free(current_solution);
	return 0;
}
