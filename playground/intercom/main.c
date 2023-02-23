#include <stdio.h>
#include <mpi.h>


#define TAG 0

int main(int argc, char* argv[])
{

  // Group1 {senders}, Group 2 {receivers}
  MPI_Comm local_comm; /* communicator for each group */
  MPI_Comm snd2rcv_comm; /* inter-communicator */
  MPI_Status status;
  int nreceivers = 2;
  int myrank, membershipKey, my_local_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


  /* recievers have membershipkey = 0 */
  membershipKey = (myrank < nreceivers);

  /* Build intra-communicator for local sub-group */
  MPI_Comm_split(MPI_COMM_WORLD, membershipKey, myrank, &local_comm);
  
  MPI_Comm_rank(local_comm, &my_local_rank);

  if (membershipKey == 0) {
    MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, 0, TAG, &snd2rcv_comm);
  } else {
    MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, nreceivers, TAG, &snd2rcv_comm);
  }

 


  // sending and receiving:

  int buf = myrank;

  if (my_local_rank == 1 && myrank >= nreceivers) {
    buf = myrank;

    MPI_Send(&buf, 1, MPI_INT, 1, TAG, snd2rcv_comm);
    printf("rank%d locally%d sent a message\n", myrank, my_local_rank);
  }

  if (my_local_rank == 1 && myrank < nreceivers) {

    MPI_Recv(&buf, 1, MPI_INT, 1, TAG, snd2rcv_comm, &status);
    printf("rank%d locally%d received a message %d\n", myrank, my_local_rank, buf);
    printf("received from %d, this is the local rank in the remote group\n", status.MPI_SOURCE);
  }


  // test gather

  int data[10] = {0};
  if (membershipKey == 0) {
    // just call gather
    MPI_Gather(&buf, 1, MPI_INT, NULL, 1, MPI_INT, 0, snd2rcv_comm) ;
  } else {
    if (my_local_rank == 0) {
      // I am root
      MPI_Gather(NULL, 1, MPI_INT, data, 1, MPI_INT, MPI_ROOT, snd2rcv_comm) ;
      puts("data=");
      for (int i = 0; i<10; ++i) 
	printf("%d, ", data[i]);
      puts("");
      
    }
    MPI_Gather(NULL, 1, MPI_INT, data, 1, MPI_INT, MPI_PROC_NULL, snd2rcv_comm) ;
  }




  MPI_Finalize();
  
  
}
