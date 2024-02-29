#include "tsp.h"

void tsp_free(struct tsp* tsp)
{
	free(tsp->xcoord);
	free(tsp->ycoord);
}
