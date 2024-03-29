#include "tsp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
	tsp->force_stop = 0;
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
		if (!strcmp(argv[i], "--random") || !strcmp(argv[i], "-r")) {
			if (modelSource != -1) {
				perror("Can't use more than 1 mode\n");
				return -1;
			}
			modelSource = 0;
		} else if (!strcmp(argv[i], "--seed") || !strcmp(argv[i], "-s")) {
			tsp->seed = atoi(argv[++i]);
			userSetSeed = 1;
		} else if (!strcmp(argv[i], "--nnodes") || !strcmp(argv[i], "-n")) {
			tsp->nnodes = atoi(argv[++i]);
			userSetNnodes = 1;
		} else if (!strcmp(argv[i], "--timelimit") || !strcmp(argv[i], "-t")) {
			tsp->timelimit_secs = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--inputfile") || !strcmp(argv[i], "-i")) {
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
			perror("Set seed and nnodes to use random generation\n");
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
	FILE* where = stderr;
	fprintf(where, "----------INSTANCE-----------\n");
	fprintf(where, "nnodes: %d\n", tsp->nnodes);

	fprintf(where, "model_source: ");
	if (tsp->model_source == 1)
		fprintf(where, "random\n");
	else if (tsp->model_source == 2)
		fprintf(where, "input_file\n");
	else
		fprintf(where, "NOT SET\n");

	fprintf(where, "seed: %d\n", tsp->seed);
	fprintf(where, "input_file: %s\n", tsp->input_file);
	fprintf(where, "edge weight type: %s\n", tsp->edge_weight_type);

	if (tsp->solution_permutation) {
		for (int i = 0; i < tsp->nnodes; i++) {
			fprintf(where, "%d->", tsp->solution_permutation[i]);
		}
		fprintf(where, "\nBEST SOLUTION: %lf\n", tsp->solution_value);
	}
	fprintf(where, "-----------------------------\n");
}

void debug_print_coords(struct tsp* tsp)
{
	FILE* where = stderr;
	fprintf(where, "------------NODES------------\n");
	for (int i = 0; i < tsp->nnodes; i++) {
		fprintf(where, "%d: (%lf, %lf)\n", i + 1, tsp->coords[i].x, tsp->coords[i].y);
	}
	fprintf(where, "-----------------------------\n");
}

int tsp_allocate_costs(struct tsp* tsp)
{
	if (tsp->nnodes <= 0)
		return -1;

	if (tsp->cost_matrix)
		free(tsp->cost_matrix);

	tsp->cost_matrix = (double*)malloc(sizeof(double) * tsp->nnodes * tsp->nnodes);
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

int nint(double x)
{
	return (int)(x + 0.5);
}

int tsp_compute_costs(struct tsp* tsp, tsp_costfunction costfunction)
{
	if (tsp->cost_matrix == NULL)
		tsp_allocate_costs(tsp);

	if (tsp->cost_matrix == NULL || tsp->coords == NULL)
		return -1;

	for (int i = 0; i < tsp->nnodes; i++) {
		for (int j = 0; j < tsp->nnodes; j++) {

			double dij = costfunction(tsp->coords[i].x, tsp->coords[j].x, tsp->coords[i].y,
						  tsp->coords[j].y);

			tsp->cost_matrix[flatten_coords(i, j, tsp->nnodes)] = dij;
		}
	}

	return 0;
}

double compute_delta(struct tsp* tsp, int* solution, int i, int j)
{
	double
	    distance_prev = tsp->cost_matrix[flatten_coords(solution[i], solution[i + 1], tsp->nnodes)] +
			    tsp->cost_matrix[flatten_coords(solution[j], solution[(j + 1) % tsp->nnodes], tsp->nnodes)];
	double distance_next = tsp->cost_matrix[flatten_coords(solution[i + 1], solution[(j + 1) % tsp->nnodes],
							       tsp->nnodes)] +
			       tsp->cost_matrix[flatten_coords(solution[i], solution[j], tsp->nnodes)];
	return distance_prev - distance_next;
}

int tsp_2opt_solution(struct tsp* tsp,
		      int* current_solution,
		      double* current_solution_value,
		      int best_i,
		      int best_j,
		      double best_delta)
{
	tsp_2opt_swap(best_i + 1, best_j, current_solution);
	*current_solution_value -= best_delta;
	if (*current_solution_value < tsp->solution_value) {
		tsp->solution_value = *current_solution_value;
		tsp_save_solution(tsp, current_solution, *current_solution_value);
		return 1;
	}
	return 0;
}

double tsp_recompute_solution_arg(struct tsp* tsp, int* solution)
{
	double current_solution = 0;
	for (int i = 0; i < tsp->nnodes - 1; i++) {
		current_solution += tsp->cost_matrix[flatten_coords(solution[i], solution[i + 1], tsp->nnodes)];
	}
	current_solution += tsp->cost_matrix[flatten_coords(solution[0], solution[tsp->nnodes - 1], tsp->nnodes)];

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

void tsp_save_solution(struct tsp* tsp, int* solution, double value)
{
	memcpy(tsp->solution_permutation, solution, sizeof(int) * tsp->nnodes);
	tsp->solution_value = value;
}

int tsp_is_solution_arg(int* solution, int nnodes)
{
	if (solution == NULL)
		return 0;
	char* found = calloc(nnodes, sizeof(char));

	int res = 1;

	for (int i = 0; i < nnodes; i++) {
		found[solution[i]]++;
	}

	for (int i = 0; i < nnodes; i++) {
		if (found[i] != 1) {
			res = 0;
			break;
		}
	}

	free(found);
	return res;
}

int tsp_is_solution(struct tsp* tsp)
{
	return tsp_is_solution_arg(tsp->solution_permutation, tsp->nnodes);
}

double tsp_2opt_findbestswap(struct tsp* tsp, int* solution, int* best_i, int* best_j)
{
	double best_delta = -10e30;
	for (int i = 0; i < tsp->nnodes - 2; i++) {
		for (int j = i + 2; j < tsp->nnodes; j++) {
			double delta = compute_delta(tsp, solution, i, j);
			if (delta > best_delta) {
				best_delta = delta;
				*best_i = i;
				*best_j = j;
			}
		}
	}
	return best_delta;
}

void tsp_2opt_swap(int left, int right, int* solution)
{
	while (left < right) {
		int temp = solution[left];
		solution[left] = solution[right];
		solution[right] = temp;
		left++;
		right--;
	}
}

double tsp_costfunction_att(double xi, double xj, double yi, double yj)
{
	double deltax = (xi - xj);
	double deltay = (yi - yj);
	double sqdist = deltax * deltax + deltay * deltay;
	double rij = sqrt(sqdist / 10.0);
	double tij = nint(rij);
	double dij;
	if (tij < rij)
		dij = tij + 1;
	else
		dij = tij;

	return dij;
}

double tsp_costfunction_euclidian(double xi, double xj, double yi, double yj)
{
	double deltax = (xi - xj);
	double deltay = (yi - yj);
	double sqdist = deltax * deltax + deltay * deltay;
	return sqrt(sqdist);
}

int tsp_shouldstop(struct tsp* tsp)
{
	time_t currenttime = clock();
	time_t starttime = tsp->start_time;

	if ((double)(currenttime - starttime) / CLOCKS_PER_SEC > (double)tsp->timelimit_secs)
		return 1;

	return 0;
}

void tsp_starttimer(struct tsp* tsp)
{
	tsp->start_time = clock();
}
