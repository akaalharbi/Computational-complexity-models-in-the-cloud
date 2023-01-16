#ifndef LONG_MESSAGE_RECEIVER
#define LONG_MESSAGE_RECEIVER

#include "numbers_shorthands.h"
#include "numbers_shorthands.h"
#include "hash.h"
#include "dict.h"
#include "config.h"
#include <mpi.h>

void receiver(int myrank, MPI_Comm mpi_communicator, int nsenders);


#endif 

