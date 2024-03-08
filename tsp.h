#ifndef TSP_H_
#define TSP_H_

#include <stdlib.h>

#define RANDOM_MAX_X 10000
#define RANDOM_MAX_Y 10000
#define EPSILON      1e-7

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

struct solver_parameters {
	int starting_node;
};
struct solver_stack;

typedef int(tsp_solver)(struct tsp* tsp,
			int* output_solution,
			double* output_value,
			struct solver_parameters solver_parameters,
			struct solver_stack* solver_stack);

struct solver_stack {
	struct solver_stack* next;
	tsp_solver* solver;
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

// SOLVERS
int tsp_solve_greedy(struct tsp* tsp,
		     int* output_solution,
		     double* output_value,
		     struct solver_parameters solver_parameters,
		     struct solver_stack* solver_stack);
int tsp_solve_allstartingnodes(struct tsp* tsp,
			       int* output_solution,
			       double* output_value,
			       struct solver_parameters solver_parameters,
			       struct solver_stack* solver_stack);
int tsp_solve_2opt(struct tsp* tsp,
		   int* output_solution,
		   double* output_value,
		   struct solver_parameters solver_parameters,
		   struct solver_stack* solver_stack);

int tsp_solve_save(struct tsp* tsp,
		   struct solver_parameters solver_parameters,
		   struct solver_stack* solver_stack);

double tsp_recompute_solution_arg(struct tsp* tsp, int* solution);
double tsp_recompute_solution(struct tsp* tsp);
int tsp_check_solution(struct tsp* tsp, double* computed);

#endif // TSP_H_
