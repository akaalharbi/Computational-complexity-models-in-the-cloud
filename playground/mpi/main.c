#include "mpi.h"
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <unistd.h>
#define TAG 0
#define SIZE 10

void print_array(int* a, int len){
  printf("{");
  for (int i = 0; i<len; i++) {
    printf("%d, ", a[i]);
  } puts("}");
}

int main(int argc, char* argv[]){

  int rank, size;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int snd_buf = rank;
  int rcv_buf[SIZE] = {0};
  sleep(rank);

  
  if (rank == 0){
    MPI_Request requests[SIZE];
    MPI_Status statuses[SIZE];
    
    int outcount = 0;
    int indices[SIZE] = {0};

    for (int i = 0; i<SIZE; i++){
      MPI_Irecv(rcv_buf+i, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[i]);
    }
    

    int ctr=0;
    int flag=0;
    while(!flag){
      printf("ctr=%d\n", ctr);
      MPI_Waitsome(SIZE, requests, &outcount, indices, statuses);

      printf("outcount=%d\n", outcount);
      puts("indices");
      print_array(indices, outcount);
      puts("recieve buffer");
      print_array(rcv_buf, SIZE);
      MPI_Testall(SIZE, requests, &flag, statuses);
      printf("flag=%d\n", flag);
      ctr++;

    }
  } else {
    snd_buf = rank+1;
    //sleep(rank%4);
    MPI_Send(&snd_buf, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
  }
  

  MPI_Finalize(); 

  return 0;
}
