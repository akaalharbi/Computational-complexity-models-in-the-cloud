#ifndef TIMING_FUNCTION
#define TIMING_FUNCTION

#include <sys/time.h>
#include <bits/types/struct_timeval.h>



double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, 0);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

#endif
