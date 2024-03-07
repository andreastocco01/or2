#include "tsp.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

void tsp_free(struct tsp* tsp)
{
	if (tsp->coords)
		free(tsp->coords);

	if (tsp->edge_weight_type)
		free(tsp->edge_weight_type);

	if (tsp->cost_matrix)
		free(tsp->cost_matrix);

	if (tsp->solution_permutation)
		free(tsp->solution_permutation);
}

void tsp_init(struct tsp* tsp)
{
	memset(tsp, 0, sizeof(struct tsp));
	tsp->input_file = NULL;
	tsp->model_source = 0;
	tsp->coords = NULL;
	tsp->edge_weight_type = NULL;
	tsp->nnodes = 0;
}

int tsp_allocate_buffers(struct tsp* tsp)
{
	if (tsp->coords)
		free(tsp->coords);

	if (tsp->nnodes <= 0)
		return -1;

	tsp->coords = (struct point*)malloc(sizeof(struct point) * tsp->nnodes);

	return 0;
}

/*
 * Parses command line arguments for the tsp instance
 * returns:
 * -1 if parse fails
 * 0 otherwise
 * */
int tsp_parse_arguments(int argc, char** argv, struct tsp* tsp)
{
	int modelSource = -1; // 0 random, 1 file load
	int userSetSeed = 0;
	int userSetNnodes = 0;
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-random")) {
			if (modelSource != -1) {
				perror("Can't use more than 1 mode\n");
				return -1;
			}
			modelSource = 0;
		} else if (!strcmp(argv[i], "-seed")) {
			tsp->seed = atoi(argv[++i]);
			userSetSeed = 1;
		} else if (!strcmp(argv[i], "-nnodes")) {
			tsp->nnodes = atoi(argv[++i]);
			userSetNnodes = 1;
		} else if (!strcmp(argv[i], "-inputfile")) {
			tsp->input_file = argv[++i];
			if (modelSource != -1) {
				perror("Can't use more than 1 mode\n");
				return -1;
			}
			modelSource = 1;
		}
	}

	if (modelSource == 0) {
		if (userSetSeed == 0 | userSetNnodes == 0) {
			perror(
			    "Set seed and nnodes to user random generation\n");
			return -1;
		}
		tsp->model_source = 1;
	}

	else if (modelSource == 1) {
		if (tsp->input_file == NULL) {
			perror("Set an input file to use model loading\n");
			return -1;
		}
		tsp->model_source = 2;
	} else if (modelSource == -1) {
		perror("Please specify a mode to load the file\n");
		return -1;
	}

	return 0;
}

void debug_print(struct tsp* tsp)
{
	printf("nnodes: %d\n", tsp->nnodes);

	printf("model_source: ");
	if (tsp->model_source == 1)
		printf("random\n");
	else if (tsp->model_source == 2)
		printf("input_file\n");
	else
		printf("NOT SET\n");

	printf("seed: %d\n", tsp->seed);
	printf("input_file: %s\n", tsp->input_file);
	printf("edge weight type: %s\n", tsp->edge_weight_type);

	if (tsp->solution_permutation) {
		for (int i = 0; i < tsp->nnodes; i++) {
			printf("%d->", tsp->solution_permutation[i]);
		}
		printf("\nBEST SOLUTION: %lf\n", tsp->solution_value);
	}
}

void debug_print_coords(struct tsp* tsp)
{
	printf("%d nodes:\n", tsp->nnodes);
	for (int i = 0; i < tsp->nnodes; i++) {
		printf("%d: (%lf, %lf)\n", i + 1, tsp->coords[i].x,
		       tsp->coords[i].y);
	}
}

int tsp_allocate_costs(struct tsp* tsp)
{
	if (tsp->nnodes <= 0)
		return -1;

	if (tsp->cost_matrix)
		free(tsp->cost_matrix);

	tsp->cost_matrix = (double*)malloc(sizeof(double) * tsp->nnodes *
					   tsp->nnodes);
	return 0;
}

int tsp_allocate_solution(struct tsp* tsp)
{
	if (tsp->nnodes <= 0)
		return -1;

	if (tsp->solution_permutation)
		free(tsp->solution_permutation);

	tsp->solution_permutation = (int*)malloc(sizeof(int) * tsp->nnodes);

	return 0;
}

#define flatten_coords(x, y, N) x* N + y
int tsp_compute_costs(struct tsp* tsp)
{
	if (tsp->cost_matrix == NULL || tsp->coords == NULL)
		return -1;

	for (int i = 0; i < tsp->nnodes; i++) {
		for (int j = 0; j < tsp->nnodes; j++) {
			double deltax = (tsp->coords[i].x - tsp->coords[j].x);
			double deltay = (tsp->coords[i].y - tsp->coords[j].y);
			double sqdist = deltax * deltax + deltay * deltay;

			tsp->cost_matrix[flatten_coords(
			    i, j, tsp->nnodes)] = sqrt(sqdist);
		}
	}

	return 0;
}

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

// TODO move to util
void print_array(int* arr, int size)
{
	for (int i = 0; i < size; i++) {
		printf("%d->", arr[i]);
	}
	printf("\n");
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
		printf("Solving %d/%d\n", i + 1, tsp->nnodes);
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

int tsp_2opt_solution(struct tsp* tsp, int* solution, double* output_value)
{
	if (!tsp->cost_matrix)
		return -1;

	if (!tsp->nnodes)
		return -1;
start:
	for (int i = 0; i < tsp->nnodes - 2; i++) {
		for (int j = i + 2; j < tsp->nnodes; j++) {
			double
			    distance_prev = tsp->cost_matrix[flatten_coords(
						solution[i], solution[i + 1],
						tsp->nnodes)] +
					    tsp->cost_matrix[flatten_coords(
						solution[j],
						solution[(j + 1) % tsp->nnodes],
						tsp->nnodes)];
			double
			    distance_next = tsp->cost_matrix[flatten_coords(
						solution[i + 1],
						solution[(j + 1) % tsp->nnodes],
						tsp->nnodes)] +
					    tsp->cost_matrix[flatten_coords(
						solution[i], solution[j],
						tsp->nnodes)];

			double delta = distance_prev - distance_next;

			if (delta > 0) {
				int left = i + 1;
				int right = j;
				while (left < right) {
					int temp = solution[left];
					solution[left] = solution[right];
					solution[right] = temp;
					left++;
					right--;
				}

				double old_value = *output_value;
				*output_value -= delta;
				/* *output_value = tsp_recompute_solution_arg(
				 */
				/*     tsp, solution); */

				/* if (fabs((old_value - *output_value) - delta)
				 * >= */
				/*     EPSILON) { */
				/* 	printf("inconsistent!!!!!!\n"); */
				/* } */
				goto start;
			}
		}
	}
	return 0;
}

double tsp_recompute_solution_arg(struct tsp* tsp, int* solution)
{
	double current_solution = 0;
	for (int i = 0; i < tsp->nnodes - 1; i++) {
		current_solution += tsp->cost_matrix[flatten_coords(
		    solution[i], solution[i + 1], tsp->nnodes)];
	}
	current_solution += tsp->cost_matrix[flatten_coords(
	    solution[0], solution[tsp->nnodes - 1], tsp->nnodes)];

	return current_solution;
}

double tsp_recompute_solution(struct tsp* tsp)
{
	return tsp_recompute_solution_arg(tsp, tsp->solution_permutation);
}

/**
 * 1 if correct
 * 0 otherwise
 *
 * if computed is set, return the correct solution
 * */
int tsp_check_solution(struct tsp* tsp, double* computed)
{
	double computedsol = tsp_recompute_solution(tsp);

	if (computed)
		*computed = computedsol;

	if (fabs(computedsol - tsp->solution_value) >= EPSILON) {
		return 0;
	}
	return 1;
}
