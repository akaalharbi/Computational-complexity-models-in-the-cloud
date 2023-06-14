#include "numbers_shorthands.h"
#include "hash.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h> // access functoin
#include "config.h"
#include "timing.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include "util_files.h"
#include "common.h"
#include "c_sha256_avx.h"
#include <mpi.h>

void write_states(u32 tr_states[restrict 16 * 8], /* assume avx512 */
                  FILE* files[16],
		  int n_activ_lanes /* we may have less than 16 active lanes */
		  )
{
  /// take a transposed state, write each state in a file
  /// file names: 16*rank + i where i=0, 1, ..., 15
  /// then a shell script to combine all states

  u32 states[NWORDS_STATE*16];
  untranspose_state(states, tr_states);
  
  for (int i=0; i<n_activ_lanes; ++i)
    /* save the ith state in the ith file */
    fwrite(&states[i*8], 32, 1, files[i]);
}



void split_states(int myrank,
		  int nprocesses,
		  MPI_Comm inter_comm)
{
  // ==========================================================================+
  // Summary: split the state to make the interval 2^23
  // --------------------------------------------------------------------------+
  FILE* fp = fopen("data/states", "r");
  const size_t new_interval = (1<<23);
  
  u8 Mavx[16][HASH_INPUT_SIZE] = {0};
  u32 tr_states[16*8] = {0}; /* same as current_states but transposed */
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = (myrank * nstates)/nprocesses;
  size_t end = ((myrank + 1) * nstates)/nprocesses;
  int n_active_lanes;
  
  /* we would like (end - begin) = 1 + 16y to avoid problem with boundary */
  end = end + ((1 + ((begin-end)%16)) % 16);
  
  size_t global_idx; /* where are we in the states file */
  size_t local_idx; /* where are we in the buffer copied from states file  */
  int inited = 0; /* 0 if we need to clear the avx register */


  /* save my output here */
  FILE* files[16];
  char file_name[FILE_NAME_MAX_LENGTH]; /* "data/messages/%d" */
  for (int i=0; i<16; ++i) {
    /* create files */
    snprintf(file_name, sizeof(file_name), "data/states%d", myrank*16+i );
    files[i] = fopen(file_name, "w");
  }



  /* get all states that I should work on: */
  WORD_TYPE* states = (WORD_TYPE*) malloc((end - begin)
					  * sizeof(WORD_TYPE)
					  * NWORDS_STATE);


  printf("rank=%d, begin=%lu, end=%lu, ndigests=%lu, quota=%lu\n",
	 myrank, begin, end, nstates, (end - begin));

  /* only load states that i am going to work on */
  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, HASH_STATE_SIZE, (end - begin), fp);
  
  // ---> 
  /* Hash the long message again, 16 at a time */
  for (global_idx = begin; global_idx < end-1; global_idx += 16){
    /* local_idx = 0 -> (end-global)/16 */
    local_idx = global_idx - begin ;
    n_active_lanes = MIN((end - global_idx), 16);
    inited = 0; /* please clear the avx register */
    
    /* form the state to be accepted to the uint32_t *sha256_multiple_x16_tr */
    transpose_state(tr_states, &states[local_idx*NWORDS_STATE]); // this the important
    /* untranspose_state(current_states, tr_states); */ // no need to it.


    /* set message counters */
    for (int lane = 0; lane<16; ++lane)
      ((u64*) Mavx[lane])[0] = INTERVAL * (global_idx + lane);


    for (size_t hash_n=0; hash_n < INTERVAL; ++hash_n){
      
      /* hash 16 messages and copy it to tr_states  */
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states, inited),
	     16*HASH_STATE_SIZE);
      inited = 1;
      

      /* hash_single(state_singe, Mavx[0]); */
      /* update message counters */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;

      if (hash_n % new_interval == 0)
	write_states(tr_states, files, n_active_lanes);
    } /* end for hash interval */
  }

  free(states);
  fclose(fp);
}







int main(int argc, char *argv[])
{
  int nproc, myrank;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  printf("going to call split states\n");
  split_states(myrank, nproc, MPI_COMM_WORLD);
  
  return 0;
}
