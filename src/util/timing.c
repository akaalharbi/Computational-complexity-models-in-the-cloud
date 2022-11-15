#include "timing.h"


double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, 0);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}


