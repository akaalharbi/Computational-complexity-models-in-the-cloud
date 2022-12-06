#include "mpi.h"
#include <stdio.h>
#include <stdint.h>
#include <omp.h>


int main(int argc, char* argv[]){

  int rank, size;

  MPI_INIT(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);



  MPI_Finalize(); 

  return 0;
}
