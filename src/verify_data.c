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
  u32 tmp[16*8] = {0}; /* are the states after each hashing */
  u32 next_states[16*8] = {0}; /* next in the sense after INTERVAL hashin */
  u32 tr_states[16*8] = {0}; /* same as current_states but transposed */


  u32 state_singe[8];


  
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = myrank * (nstates/nprocesses);
  size_t end = (myrank + 1) * (nstates/nprocesses);
  size_t global_idx; /* where are we in the states file */
  size_t local_idx; /* where are we in the buffer copied from states file  */
  
  int print_rank = 0; /* only for debugging */
    
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

  // let's say everything up to this point is perfect!
  printf("rank=%d, begin=%lu\n", myrank, begin);


  /* Hash the long message again, 16 at a time */
  for (global_idx = begin; global_idx < end; global_idx += 16){
    /* local_idx = 0 -> (end-global)/16 */
    local_idx = global_idx - begin ;

    /* form the state to be accepted to the uint32_t *sha256_multiple_x16_tr(u32*, u32*) */
    transpose_state(tr_states, &states[local_idx*NWORDS_STATE]);

    
    print_byte_txt("tr_state", (u8*)tr_states, HASH_STATE_SIZE*16);
    puts("");    
    untranspose_state(current_states, tr_states);
    print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16);
    puts("");
    print_byte_txt("state", (u8*)&states[local_idx*16], HASH_STATE_SIZE*16);



    
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
      /* hash 16 messages and copy it to tr_states  */
      memcpy(tmp,
	     sha256_multiple_x16_tr(Mavx, tr_states),
	     16*HASH_STATE_SIZE);
      
      memcpy(tr_states, tmp, 16*HASH_STATE_SIZE);

      /* hash_single(state_singe, Mavx[0]); */
      /* update message counters */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;


      /* /\* check we have the same hashes *\/ */
      /* untranspose_state(current_states, tr_states); */

      /* nbytes_non_equal = memcmp(current_states, state_singe, HASH_STATE_SIZE); */
      /* if (nbytes_non_equal != 0) { */
      /* 	printf("hurray at %lu\n", hash_n); */
      /* 	print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16); */
      /* 	puts(""); */
      /* 	print_byte_txt("single", (u8*)state_singe, HASH_STATE_SIZE); */
      /* 	puts(""); */
      /* } */

      
    } /* end for hash interval */


    /* check we have the same hashes */
    untranspose_state(current_states, tr_states);

    nbytes_non_equal = memcmp(current_states, next_states, 16*HASH_STATE_SIZE);
    is_corrupt = (0 != nbytes_non_equal);

    if (is_corrupt) {
      printf("found a curropt state at global_idx=%lu\n", global_idx);
      printf("first hash=%d\n",
	     memcmp(current_states, next_states, 16*HASH_STATE_SIZE));
      
      fprintf(fp_verify, "found a curropt state at global_idx=%lu\n", global_idx);
      print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16);
      puts("");
      print_byte_txt("next", (u8*)next_states, HASH_STATE_SIZE*16);
    }
    elapsed = wtime() - elapsed;
    printf("rank=%d, step=%lu, is_corrupt=%d, elapsed=%0.2fsec\n",
	   myrank,
	   global_idx,
	   is_corrupt,
	   elapsed);

    
    return; /* just one hash */
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
