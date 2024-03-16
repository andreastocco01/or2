#ifndef TSP_H_
#define TSP_H_

#define RANDOM_MAX_X            10000
#define RANDOM_MAX_Y            10000
#define EPSILON                 1e-7
#define flatten_coords(x, y, N) x* N + y

#include <time.h>

struct point {
	double x;
	double y;
};

struct tsp {
	// instance data
	int nnodes;
	struct point* coords;

	int model_source; // 1 if random, 2 if input_file

	// random generator data
	int seed;

	// input file data
	char* input_file;

	char* edge_weight_type;

	double* cost_matrix;

	int* solution_permutation;
	double solution_value;
};

void tsp_init(struct tsp* tsp);
int tsp_allocate_buffers(struct tsp* tsp);
int tsp_allocate_solution(struct tsp* tsp);
int tsp_allocate_costs(struct tsp* tsp);
void tsp_free(struct tsp* tsp);
int tsp_parse_arguments(int argc, char** argv, struct tsp* tsp);
void debug_print(struct tsp* tsp);
void debug_print_coords(struct tsp* tsp);

int tsp_compute_costs(struct tsp* tsp);

int tsp_2opt_solution(struct tsp* tsp, int* solution, double* output_value);

double tsp_recompute_solution_arg(struct tsp* tsp, int* solution);
double tsp_recompute_solution(struct tsp* tsp);
int tsp_check_solution(struct tsp* tsp, double* computed);

void tsp_save_signal_safe(struct tsp* tsp, int* solution, double value);

#endif // TSP_H_
