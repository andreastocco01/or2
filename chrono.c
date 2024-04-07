#include "chrono.h"
#include <time.h>
#include <stdlib.h>

#include <sys/resource.h>

double walltime()
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + 1.0e-9*((double)ts.tv_nsec);
}

double second()
{
	double t = walltime();
	return(t);
}
