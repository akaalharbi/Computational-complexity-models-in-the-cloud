#ifndef LONG_MESSAGE_SENDER
#define LONG_MESSAGE_SENDER


#include "numbers_shorthands.h"
#include "config.h"
#include <mpi.h>


void sender(int myrank, MPI_Comm mpi_communicator);


#endif
