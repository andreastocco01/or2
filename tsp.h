#ifndef TSP_H_
#define TSP_H_

#include <math.h>

#define RANDOM_MAX_X               10000
#define RANDOM_MAX_Y               10000
#define EPSILON                    1e-7
#define STARTING_INCUMBENTS        4096
#define STARTING_CURRENT_SOLUTIONS 4096
#define flatten_coords(x, y, N)    x* N + y

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

// COST FUNCTIONS
typedef double (*tsp_costfunction)(double xi, double xj, double yi, double yj);

int nint(double x);

double tsp_costfunction_att(double xi, double xj, double yi, double yj);
double tsp_costfunction_euclidian(double xi, double xj, double yi, double yj);

/**
 * Initialize a tsp structure
 * */
void tsp_init(struct tsp* tsp);

/**
 * Allocate buffers of a tsp structure.
 *
 * IMPORTANT: nnodes must be set
 * */
int tsp_allocate_buffers(struct tsp* tsp);

/**
 * Allocate the data structure used to store the solution
 * */
int tsp_allocate_solution(struct tsp* tsp);

/**
 * Allocate the data structure used to save the costs
 * */
int tsp_allocate_costs(struct tsp* tsp);

/**
 * Fills the matrix of costs
 * */
int tsp_compute_costs(struct tsp* tsp, tsp_costfunction costfunction);

/**
 * Free memory allocated by a tsp struct
 * */
void tsp_free(struct tsp* tsp);

/**
 * Parse the standard argument passed as line arguments
 * */
int tsp_parse_arguments(int argc, char** argv, struct tsp* tsp);

void debug_print(struct tsp* tsp);
void debug_print_coords(struct tsp* tsp);

/**
 * Compute the delta that would be obtained by applying
 * the 2opt procedure on the given indices
 * */
double compute_delta(struct tsp* tsp, int* solution, int i, int j);

/**
 * Executes a 2-opt swap and updates the current solution.
 *
 * returns 1 if a new best solution (incumbent) is found
 * */
int tsp_2opt_solution(struct tsp* tsp,
		      int* current_solution,
		      double* current_solution_value,
		      int best_i,
		      int best_j,
		      double best_delta);

/**
 * Finds the best 2opt swap nodes
 *
 * returns the value of the delta
 * */
double tsp_2opt_findbestswap(struct tsp* tsp, int* solution, int* best_i, int* best_j);

/**
 * Execute the 2opt swap operation on a solution
 * */
void tsp_2opt_swap(int left, int right, int* solution);

/**
 * Compute the cost of a solution passed as argument
 * */
double tsp_recompute_solution_arg(struct tsp* tsp, int* solution);

/**
 * Compute the cost of the solution inside the struct
 * */
double tsp_recompute_solution(struct tsp* tsp);

/**
 * Compute the cost of the solution again and check if it corresponds
 * with the one saved inside the struct.
 *
 * computed: output of the value
 * */
int tsp_check_solution(struct tsp* tsp, double* computed);

/**
 * Check whether the solution passed as argument is a valid cycle
 *
 * returns 1 if it is, 0 otherwise
 * */
int tsp_is_solution_arg(int* solution, int nnodes);

/**
 * Check whether the solution in the tsp struct is valid cycle
 *
 * returns 1 if it is, 0 otherwise
 * */
int tsp_is_solution(struct tsp* tsp);

/**
 * Save a new solution inside the tsp struct in a signal-safe way
 * */
void tsp_save_signal_safe(struct tsp* tsp, int* solution, double value);

#endif // TSP_H_
