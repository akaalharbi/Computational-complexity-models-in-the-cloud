#ifndef LONG_MESSAGE_SENDER
#define LONG_MESSAGE_SENDER


#include "numbers_shorthands.h"
#include "config.h"
#include <mpi.h>



void sender(MPI_Comm local_comm, MPI_Comm inter_comm);

#endif
