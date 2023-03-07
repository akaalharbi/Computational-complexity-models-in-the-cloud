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

void verify_region(size_t start, size_t end)
{

  #pragma omp parallel for
  for (size_t state_number=start; state_number<end; ++state_number) {
    u64 state_number = 5101;

    u64 ctr = state_number*INTERVAL;
  
    WORD_TYPE state_before[NWORDS_STATE];
    WORD_TYPE state_after[NWORDS_STATE];
  
    FILE* fp = fopen("data/states", "r");
    size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
    printf("There are %lu middle states\n", nstates);
  
    fseek(fp, (HASH_STATE_SIZE)*(state_number), SEEK_CUR);

  
    fread(state_before, 1, HASH_STATE_SIZE,  fp);
    fread(state_after, 1,  HASH_STATE_SIZE, fp);
    fclose(fp);
    
    print_char((u8*) state_before, HASH_STATE_SIZE);
    print_char((u8*) state_after, HASH_STATE_SIZE);
    printf("ctr=%llu, INTERVAL=%llu\n", ctr, INTERVAL);
    puts("Going to check...");


  

    int nbytes_non_equal = check_hashes_interval_single(state_before,
							state_after,
							ctr);

    int is_corrupt = (0 != nbytes_non_equal);

    printf("nbytes_non_equal=%d\n", nbytes_non_equal);
    printf("Found a corrupt state? %d\n", is_corrupt);

  }

  
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
  /* u32 tmp[16*8] = {0}; /\* are the states after each hashing *\/ */
  u32 next_states[16*8] = {0}; /* next in the sense after INTERVAL hashin */
  u32 tr_states[16*8] = {0}; /* same as current_states but transposed */


  u32 state_singe[8];


  
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = myrank * (nstates/nprocesses);
  size_t end = (myrank + 1) * (nstates/nprocesses);
  size_t global_idx; /* where are we in the states file */
  size_t local_idx; /* where are we in the buffer copied from states file  */
  int inited = 0; /* 0 if we need to clear the avx register */
 
    
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


  /* Hash the long message again, 16 at a time */
  for (global_idx = begin; global_idx < end; global_idx += 16){
    /* local_idx = 0 -> (end-global)/16 */
    local_idx = global_idx - begin ;
    inited = 0; /* please clear the avx register */
    
    /* form the state to be accepted to the uint32_t *sha256_multiple_x16_tr */
    transpose_state(tr_states, &states[local_idx*NWORDS_STATE]);

    /* this line should be removed as well */
    untranspose_state(current_states, tr_states);
    
    memcpy(state_singe, &states[local_idx*NWORDS_STATE], HASH_STATE_SIZE);

    /* we will test eventually transpose(tr_states) =?= next_states */
    memcpy(next_states,
	   &states[(local_idx + 1)*NWORDS_STATE],
	   HASH_STATE_SIZE*16);

    /* set message counters */
    for (int lane = 0; lane<16; ++lane)
      ((u64*) Mavx[lane])[0] = INTERVAL * (global_idx + lane);


    elapsed = wtime();

    for (size_t hash_n=0; hash_n < INTERVAL; ++hash_n){
      /* we can get rid of tmp copying */
      
      /* hash 16 messages and copy it to tr_states  */
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states, inited),
	     16*HASH_STATE_SIZE);
      inited = 1;
      

      /* hash_single(state_singe, Mavx[0]); */
      /* update message counters */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;

      
    } /* end for hash interval */


    /* check we have the same hashes */
    untranspose_state(current_states, tr_states);

    nbytes_non_equal = memcmp(current_states, next_states, 16*HASH_STATE_SIZE);
    is_corrupt = (0 != nbytes_non_equal);

    if (is_corrupt && myrank==0) {
      printf("found a curropt state at global_idx=%lu\n", global_idx);
      printf("first hash=%d\n",
	     memcmp(current_states, next_states, 16*HASH_STATE_SIZE));
      
      fprintf(fp_verify, "found a curropt state at global_idx=%lu\n", global_idx);
      print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16);
      puts("");
      print_byte_txt("next", (u8*)next_states, HASH_STATE_SIZE*16);
    }
    elapsed = wtime() - elapsed;
    printf("rank=%d, step=%lu, corrupt?=%d, elapsed=%0.2fsec, 2^%0.2f\n",
	   myrank,
	   global_idx,
	   is_corrupt,
	   elapsed,
	   log2(INTERVAL/elapsed)+log2(16));
    /* 16*/
    
    /* return; /\* just one hash *\/ */
  }
  
  free(states);
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
  fclose(fp);
  printf("There are %lu middle states\n", nstates);

  puts("Going to check...");
  verify_middle_states(rank, size, MPI_COMM_WORLD);  



  MPI_Finalize();
  return 0;
  
}
