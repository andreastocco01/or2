#include "util.h"
#include <stdio.h>
#include <stdlib.h>

double random01()
{
	return ((double)rand() / RAND_MAX);
}

void print_array_int(int* arr, int size)
{
	for (int i = 0; i < size; i++) {
		printf("%d->", arr[i]);
	}
	printf("\n");
}

void print_array_double(double* arr, int size)
{
	for (int i = 0; i < size; i++) {
		printf("%lf->", arr[i]);
	}
	printf("\n");
}
