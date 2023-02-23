// phase ii: high level overview 
// four types of processors: senders (the majority), receivers (#NSERVERS),
// Upon starting the program:
//
// `receiver`:- (init) load digests into a dictionary(once), receive messages
//              templates from senders (once).
//            - (listen)recieve digests from senders, probe them, and store the
//               candidates.
//            - (send to archive) if a receiver had enough digests, send them to
//              to the archive. Then go back to the listen state.
//
// `sender`:- (init) generate random message templates, send it to all recivers
//          - (gen) hash many random messages, decides which receiver x should
//              get a sepcific digest, store the digest in a buffer x, if it
//              is reach the quota, send that buffer to receiver x.
//              repeat infinitely
//





#include "numbers_shorthands.h"
#include "hash.h"

#include "dict.h"
#include <bits/types/struct_timeval.h>

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>
#include "config.h"
#include "timing.h"
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include "util_files.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include <mpi.h>
#include "common.h"
#include "sender.h"
#include "receiver.h"

// checklist: msg||dgst 





int main(int argc, char* argv[])
{

  // the program is a bit stupid, it can't handle when N <= DEFINED_BYTES 
  assert(N - DEFINED_BYTES > 0);

  // some extra arguments to communicate with other
  
  // ==========================================================================+
  // Summary: Three main characters:                                           |
  //  1- producers: rank [NSERVERS+1, INF] generate hashes indifinitely.       |
  //  2- consumers: ranks [0, NSERVERS-1] receive hashes and probe them in the |
  //                in its dictionary. Save those that return positive answer  |
  // ==========================================================================+
  // Note: recievers' global and local rank are the same!
  //       this is not the case for the senders.


  // --------------------- INIT MPI & Shared Variables ------------------------+
  int nproc, myrank;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);


  /* Variables for the inter-communicator */
  /* color 0 group -> color 1 group i.e. */
  /* senders group send to receivers group :) */
  MPI_Comm local_comm, inter_comm;

  /* color: 1 for recievers, 0 for senders */
  int color = (myrank < NSERVERS);

  /* Create a local communicator: */
  /* (senders has a local comm), (receivers has a local comme) */
  MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &local_comm);
  
  /* How many procs that are going t send */





  // Who am I? a sender,  or a receiver?
  if (myrank >= NSERVERS){
    /* Creat inter-comm from sender point of view:  */
    /* local leader: local rank 0, remote leader: global rank NSERVES, tag=0  */
    MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, NSERVERS, 0, &inter_comm);

    /* It knows the number of receivers from NSERVERS from config.h */
    /* sender(inter_comm);  */
  }
  else if (myrank < NSERVERS){ /* receiver, repeat infinitely  */
    /* Creat inter-comm from sender point of view:  */
    /* local leader: local rank 0(global rank NSERVERS), remote leader: global rank 0, tag=0  */
    MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, 0, 0, &inter_comm);

    int nsenders, local_rank /* local rank */;
    MPI_Comm_remote_size(inter_comm, &nsenders);
    MPI_Comm_rank(local_comm, &local_rank);
    
    receiver(local_rank, nsenders, inter_comm);
  }


  // The end that will never be reached in any case!
  printf("Process #%d reached the end (should be a receiver)\n", myrank);
  MPI_Finalize();
  return 0;
}







