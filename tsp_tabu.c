#include "tsp_tabu.h"

#include "eventlog.h"
#include "tsp.h"
#include "tsp_greedy.h"
#include <bits/types/sigset_t.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int max(int a, int b)
{
	return a < b ? b : a;
}

int tenure_min = 10;

void tenure_setmin(int min)
{
	tenure_min = min;
}

int tenure_fixed_divisor = 10;

int tenure_fixed(int nnodes, int iteration)
{
	int computed = nnodes / tenure_fixed_divisor;
	return max(computed, tenure_min);
}

void tenure_fixed_setdivisor(int divisor)
{
	tenure_fixed_divisor = divisor;
}

int tenure_sin_divisor = 2000;
int tenure_sin_scale = 20;

int tenure_sin(int nnodes, int iteration)
{
	double scale = (double)nnodes / tenure_sin_scale;
	int computed = sin(((double)iteration) / tenure_sin_divisor) * scale + scale * 2;
	return max(computed, tenure_min);
}

void tenure_sin_setscale(int scale)
{
	tenure_sin_scale = scale;
}

void tenure_sin_setdivisor(int divisor)
{
	tenure_sin_divisor = divisor;
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
	eventlog_logdouble("new_current", 0, current_solution_value);
	eventlog_logdouble("new_incumbent", 0, current_solution_value);

	int current_iteration = 0;
	int ten;

	tsp_starttimer(tsp);

	while (1) {
		// Intensification phase
		int best_i = -1, best_j = -1;
		double best_delta;
		while (1) {
			if (tsp_shouldstop(tsp))
				goto free_solution_buffers;
			current_iteration++;
			ten = tenure(tsp->nnodes, current_iteration);
			best_delta = tsp_2opt_findbestswap_no_tabu(tsp, current_solution, &best_i, &best_j,
								   tabu_iteration, ten, current_iteration);
			if (best_delta <= 0)
				break; // local minimum

			int isnewbest = tsp_2opt_swap_save(tsp, current_solution, &current_solution_value, best_i,
							   best_j, best_delta);
			if (isnewbest)
				eventlog_logdouble("new_incumbent", current_iteration, current_solution_value);
			eventlog_logdouble("new_current", current_iteration, current_solution_value);

			// output this even tho it is not changed in order
			// to plot a better chart
			ten = tenure(tsp->nnodes, current_iteration);
			eventlog_logdouble("tenure", current_iteration, ten);
		}

		// Diversification phase
		ten = tenure(tsp->nnodes, current_iteration);
		eventlog_logdouble("tenure", current_iteration, ten);
		/* printf("tenure(%d) = %d\n", current_iteration, ten); */
		if (is_tabu(tabu_iteration, best_i, current_iteration, ten)) {
			// the node is tabu. We need to skip
			continue;
		}

		tsp_2opt_swap(best_i + 1, best_j, current_solution);
		current_solution_value -= best_delta;
		// add to tabu list
		tabu_iteration[best_i] = current_iteration;
		eventlog_logdouble("new_current", current_iteration, current_solution_value);
	}

free_solution_buffers:

	free(tabu_iteration);
	free(current_solution);
	return 0;
}
