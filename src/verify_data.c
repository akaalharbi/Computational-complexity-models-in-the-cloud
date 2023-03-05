#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "config.h"
#include "hash.h"
#include "numbers_shorthands.h"
#include "util_char_arrays.h"
#include "util_files.h"
#include "timing.h"
#include "arch_avx512_type1.h"
#include "common.h"
#include "c_sha256_avx.h"


/// Verify that the states file is not corrupted in parallel (MPI)

int check_hashes_interval_single(const WORD_TYPE state_befoe[NWORDS_STATE],
				 const WORD_TYPE state_after[NWORDS_STATE],
				 const CTR_TYPE ctr){


  u8 M[HASH_INPUT_SIZE] = {0};
  /* We don't want to modify arguments */
  WORD_TYPE state[NWORDS_STATE];
  memcpy(state, state_befoe, HASH_STATE_SIZE);
  
  /* register counter in M */
  ((u64*)M)[0] = ctr;

  for (size_t i=0; i<INTERVAL; ++i) {

    hash_single(state, M);
    ++( ((u64*)M)[0]) ;

    if (memcmp(state, state_after, HASH_STATE_SIZE) == 0){
      printf("i=%lu\n",i );
      print_char((u8*) state, HASH_STATE_SIZE);
      print_char((u8*) state_after, HASH_STATE_SIZE);

    }
  }

  
  return memcmp(state, state_after, HASH_STATE_SIZE);
  
}

static void verify_middle_states(int myrank,
				 int nprocesses,
				 MPI_Comm inter_comm)
{
  // ==========================================================================+
  // Summary: Check the middle states found in the file data/states.           |
  // --------------------------------------------------------------------------+
  FILE* fp = fopen("data/states", "r");
  int is_corrupt = 0;
  int nbytes_non_equal = 0;
  double elapsed = 0;
  
  u8 Mavx[16][HASH_INPUT_SIZE] = {0};
  u32 current_states[16*8] = {0}; /* are the states after each hashing */
  u32 next_states[16*8] = {0}; /* next in the sense after INTERVAL hashin */
  u32 tr_states[16*8] = {0}; /* same as current_states but transposed */


  /* u32 state_singe[8]; */


  
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = myrank * (nstates/nprocesses);
  size_t end = (myrank + 1) * (nstates/nprocesses);
  size_t global_idx; /* where are we in the states file */
  size_t local_idx; /* where are we in the buffer copied from states file  */
  
    
  if (myrank == (nprocesses-1))
    return; /* coward solution */
    /* end = nstates; */


  /* save my output here */
  char file_name[FILE_NAME_MAX_LENGTH]; /* "data/messages/%d" */
  snprintf(file_name, sizeof(file_name), "data/verify/%d", myrank );
  FILE* fp_verify = fopen(file_name, "w"); /* register message candidates here */

  /* get all states that I should work on: */
  WORD_TYPE* states = (WORD_TYPE*) malloc((end - begin)
					  * sizeof(WORD_TYPE)
					  * NWORDS_STATE);


  printf("rank=%d, begin=%lu, end=%lu, ndigests=%lu, quota=%lu\n",
	 myrank, begin, end, nstates, (end - begin));

  /* only load states that i am going to work on */
  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, HASH_STATE_SIZE, (end - begin), fp);
  fclose(fp);
  
  // let's say everything up to this point is perfect!
  printf("rank=%d, begin=%lu\n", myrank, begin);


  /* Hash the long message again, 16 at a time */
  for (global_idx = begin;
       global_idx < end;
       global_idx += NHASH_LANES)
    {
    /* local_idx = 0 -> (end-global)/16 */
    local_idx = global_idx - begin ;

    /* form the state to be accepted to the uint32_t *sha256_multiple_x16_tr(u32*, u32*) */
    transpose_state(tr_states, &states[local_idx*NWORDS_STATE]);
    /* memcpy(state_singe, &states[local_idx*NWORDS_STATE], HASH_STATE_SIZE); */

    /* we will test eventually transpose(tr_states) =?= next_states */
  
    memcpy(next_states,
	   &states[(local_idx + 1)*NWORDS_STATE],
	   HASH_STATE_SIZE*16);
    

    
    /* set message counters */
    for (int lane = 0; lane<NHASH_LANES; ++lane)
      ((u64*) Mavx[lane])[0] = INTERVAL * (global_idx + lane);

    elapsed = wtime();

    for (size_t hash_n=0; hash_n < INTERVAL; ++hash_n){
      /* hash 16 messages and copy it to tr_states  */
      #ifdef  __AVX512F__
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states),
	     16*HASH_STATE_SIZE);
      #endif

      #ifndef  __AVX512F__
      #ifdef    __AVX2__
      /* sha256_multiple_oct_tr(Mavx, tr_states); */

      memcpy(tr_states,
	     sha256_multiple_oct_tr(Mavx, tr_states),
	     8*HASH_STATE_SIZE);
      #endif
      #endif


      
      /* hash_single(state_singe, Mavx[0]); */


    /* /\* Hash multiple messages at once *\/ */

    /* /\* HASH 16 MESSAGES AT ONCE *\/ */
    /* tr_states = sha256_multiple_x16(Mavx);   */
    /* #endif */

    /* #ifndef  __AVX512F__ */
    /* #ifdef    __AVX2__ */
    /* /\* HASH 16 MESSAGES AT ONCE *\/ */
    /* tr_states = sha256_multiple_oct(Mavx); */
    /* #endif */
    /* #endif */
      
      /* update message counters */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;


      /* /\* check we have the same hashes *\/ */
      /* untranspose_state(current_states, tr_states); */

      /* nbytes_non_equal = memcmp(current_states, state_singe, HASH_STATE_SIZE); */
      /* if (nbytes_non_equal != 0) { */
      /* 	printf("hurray at %lu\n", hash_n); */
      /*  	print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16); */
      /* 	puts(""); */
      /* 	print_byte_txt("single", (u8*)state_singe, HASH_STATE_SIZE); */
      /* 	puts(""); */
      /* } */
    } /* end for hash interval */


    /* check we have the same hashes */
    untranspose_state(current_states, tr_states);

    nbytes_non_equal = memcmp(current_states,
			      next_states,
			      NHASH_LANES*HASH_STATE_SIZE);
    is_corrupt = (0 != nbytes_non_equal);

    /* if (is_corrupt) { */
    /*   printf("found a curropt state at global_idx=%lu\n", global_idx); */
    /*   printf("first hash=%d\n", */
    /* 	     memcmp(current_states, next_states, 16*HASH_STATE_SIZE)); */
      
    /*   fprintf(fp_verify, "found a curropt state at global_idx=%lu\n", global_idx); */
    /*   fflush(fp_verify); */
    /*   print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16); */
    /*   puts(""); */
    /*   print_byte_txt("next", (u8*)next_states, HASH_STATE_SIZE*16); */
    /* } */
    elapsed = wtime() - elapsed;
    
    printf("rank=%d, step=%lu, is_corrupt=%d, elapsed=%0.2fsec, 2^%0.2f hashes/sec\n",
	    myrank,
	    global_idx,
	    is_corrupt,
	    elapsed,
	    log2(INTERVAL/elapsed));

    fprintf(fp_verify, "rank=%d, step=%lu, is_corrupt=%d, elapsed=%0.2fsec, 2^%0.2f hashes/sec\n",
	   myrank,
	   global_idx,
	   is_corrupt,
	    elapsed,
	    log2(INTERVAL/elapsed));
    fflush(fp_verify);
  }

  fclose(fp);
  fclose(fp_verify);
}





int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  int rank, size;  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE* fp = fopen("data/states", "r");
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  printf("There are %lu middle states\n", nstates);

  puts("Going to check...");
  verify_middle_states(rank, size, MPI_COMM_WORLD);  

  fclose(fp);

  MPI_Finalize();
  return 0;
  
}
