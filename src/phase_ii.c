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



  // --------------------- INIT MPI & Shared Variables ------------------------+
  int nproc, myrank;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* How many procs that are going t send */
  int nsenders = nproc - NSERVERS;
  /* int nreceivers = NSERVERS; */
  printf("There are %d senders\n", nsenders);


  // Who am I? a sender,  or a receiver?
  if (myrank >= NSERVERS){
    sender(myrank, MPI_COMM_WORLD); /* never ends :) */
  }
  else if (myrank < NSERVERS){ /* receiver, repeat infinitely  */
    receiver(myrank, MPI_COMM_WORLD, nsenders);
  }


  // The end that will never be reached in any case!
  printf("Process #%d reached the end (should be a receiver)\n", myrank);
  MPI_Finalize();
  return 0;
}







