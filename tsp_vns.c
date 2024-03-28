#include "tsp_vns.h"
#include "eventlog.h"
#include "tsp.h"
#include "tsp_greedy.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UPPER 10
#define LOWER 2

int compar(const void* a, const void* b)
{
	return *((int*)a) - *((int*)b);
}

void generate_3opt_positions(struct tsp* tsp, int* positions)
{
	int size = 3;
	while (1) {
		for (int i = 0; i < size; i++)
			positions[i] = rand() % tsp->nnodes;

		qsort(positions, size, sizeof(int), compar);

		if (positions[1] > positions[0] + 1 && positions[2] > positions[1] + 1 &&
		    positions[2] + 1 < tsp->nnodes) {
			break;
		}
	}
}

void tsp_3opt_swap(int* positions, int* current_solution, int* new_solution, int size)
{
	memcpy(new_solution, current_solution, sizeof(int) * (positions[0] + 1)); // elements from 0 to i

	int pos = positions[0] + 1;
	int z = positions[1];
	// elements from j to i + 1
	while (z >= positions[0] + 1) {
		new_solution[pos++] = current_solution[z--];
	}

	z = positions[2];
	while (z >= positions[1] + 1) {
		new_solution[pos++] = current_solution[z--];
	}

	memcpy(new_solution + pos, current_solution + positions[2] + 1, sizeof(int) * (size - positions[2] - 1));
}

double compute_solution_value(struct tsp* tsp, int* solution)
{
	double solution_value = 0;
	for (int i = 0; i < tsp->nnodes - 1; i++) {
		solution_value += tsp->cost_matrix[flatten_coords(solution[i], solution[i + 1], tsp->nnodes)];
	}
	solution_value += tsp->cost_matrix[flatten_coords(solution[0], solution[tsp->nnodes - 1], tsp->nnodes)];
	return solution_value;
}

int tsp_solve_vns(struct tsp* tsp)
{
	if (tsp_allocate_solution(tsp))
		return -1;

	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;

	// starting solution is solved with a greedy approach
	int starting_node = rand() % tsp->nnodes;

	fprintf(stderr, "starting from node %d\n", starting_node);
	tsp_solve_greedy(tsp, starting_node, tsp->solution_permutation, &tsp->solution_value);

	int* current_solution = malloc(tsp->nnodes * sizeof(int));
	double current_solution_value = tsp->solution_value;

	memcpy(current_solution, tsp->solution_permutation, sizeof(int) * tsp->nnodes);
	eventlog_logdouble("new_current", 0, current_solution_value);
	eventlog_logdouble("new_incumbent", 0, current_solution_value);

	int* new_solution = (int*)malloc(tsp->nnodes * sizeof(int));

	int current_iteration = 0;

	while (1) {
		// Intensification phase
		while (1) {
			current_iteration++;
			int best_i, best_j;
			double best_delta = tsp_2opt_findbestswap(tsp, current_solution, &best_i, &best_j);

			if (best_delta <= 0)
				break; // local minimum

			int isnewbest = tsp_2opt_solution(tsp, current_solution, &current_solution_value, best_i, best_j,
							best_delta);
			if (isnewbest)
				eventlog_logdouble("new_incumbent", current_iteration, current_solution_value);
			eventlog_logdouble("new_current", current_iteration, current_solution_value);
		}

		// Diversification phase
		int positions[3];
		int r = (rand() % (UPPER - LOWER + 1)) + LOWER;
		// do 3opt [UPPER, LOWER] times
		for (int i = 0; i < r; i++) {
			generate_3opt_positions(tsp, positions); // i -> j -> k
			tsp_3opt_swap(positions, current_solution, new_solution, tsp->nnodes);
			memcpy(current_solution, new_solution, sizeof(int) * tsp->nnodes);
		}
		// now we are considering all the kicks as a single move.
		// consider moving the instruction below into the cycle to count them as
		// different moves
		current_solution_value = compute_solution_value(tsp, current_solution);
		eventlog_logdouble("new_current", current_iteration, current_solution_value);
	}
	free(new_solution);
	free(current_solution);
	return 0;
}
