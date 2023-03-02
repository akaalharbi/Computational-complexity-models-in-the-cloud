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
  u32 current_states[16*8]={0}, next_states[16*8]={0}, tmp[16*8];


  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = myrank * (nstates/nprocesses);
  size_t end = (myrank + 1) * (nstates/nprocesses);

  
  
  if (myrank == (nprocesses-1))
    return; /* coward solution */
    /* end = nstates; */


  char file_name[FILE_NAME_MAX_LENGTH]; /* "data/messages/%d" */
  snprintf(file_name, sizeof(file_name), "data/verify/%d", myrank );
  FILE* fp_verify = fopen(file_name, "w"); /* register message candidates here */

  /* get all states that i should work on: */
  WORD_TYPE* states = (WORD_TYPE*) malloc((end - begin)
					  * sizeof(WORD_TYPE)
					  * NWORDS_STATE);
  
  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, HASH_STATE_SIZE, (end - begin), fp);

  
  /* Hash buffers init */
  SHA256_ARGS args;

  int print_rank = 1;


  printf("rank=%d, begin=%lu\n", myrank, begin);
  /* Hash the long message again, 16 at a time */
  for (size_t step = 0; step<(end - begin)/16; ++step){
    memcpy(next_states,
	   &states[NWORDS_STATE*(step*16+1)],
	   HASH_STATE_SIZE*16);

    /* Read fresh states, and put them in transposed way  */
    elapsed = wtime();
    
    for (int lane = 0; lane < 16; lane++) {
      args.digest[lane + 0*16] = states[NWORDS_STATE*(step*16 + lane) + 0];
      args.digest[lane + 1*16] = states[NWORDS_STATE*(step*16 + lane) + 1];
      args.digest[lane + 2*16] = states[NWORDS_STATE*(step*16 + lane) + 2];
      args.digest[lane + 3*16] = states[NWORDS_STATE*(step*16 + lane) + 3];
      args.digest[lane + 4*16] = states[NWORDS_STATE*(step*16 + lane) + 4];
      args.digest[lane + 5*16] = states[NWORDS_STATE*(step*16 + lane) + 5];
      args.digest[lane + 6*16] = states[NWORDS_STATE*(step*16 + lane) + 6];
      args.digest[lane + 7*16] = states[NWORDS_STATE*(step*16 + lane) + 7];


      /* adjust the counter part */
      /* if (myrank == print_rank && lane == 0) { */
      /* 	printf("offset = %lu\n", NWORDS_STATE*(step*16 + lane)); */
      /* printf("BEFORE rank=%d, args.data[0] = %llu, Mavx[0] = %llu\n", */
      /* 	     myrank, ((u64*) args.data_ptr[0])[0], */
      /* 	     ((u64*) Mavx[0])[0]); */
      /* } */

      ((u64*) Mavx[lane])[0] = (begin + step*16 + lane)*INTERVAL;
      args.data_ptr[lane] = Mavx[lane]; /* point data */

      /* if (myrank == print_rank && lane == 0) { */
      /* printf("AFTER rank=%d, args.data[0] = %llu, Mavx[0] = %llu\n", */
      /* 	     myrank, ((u64*) args.data_ptr)[0], */
      /* 	     	     ((u64*) Mavx[0])[0]); */

      /* } */

    }


    for (size_t hash_n=0; hash_n < INTERVAL; ++hash_n){
      /* hash 16 messages  */
      call_sha256_x16_avx512_from_c(&args, 1);
      /* update message counters */
      for (int lane = 0; lane<16; ++lane){
	((u64*) Mavx[lane])[0] += 1;
	args.data_ptr[lane] = Mavx[lane]; /* not sure this is necessary */
      } /* end for lane  */

      /* if (myrank == pint_rank){ */
	
      /* copy_transposed_state(current_states, states, 0); */
      /* print_char((u8*) current_states, HASH_STATE_SIZE); */
      /* printf("rank=%d, args.data[0] = %llu, Mavx[0] = %llu\n", */
      /* 	     myrank, ((u64*) args.data_ptr)[0], */
      /* 	     ((u64*) Mavx)[0]); */

      /* } */

    } /* end for hash interval */


    /* check we have the same hashes */
    for (int i = 0; i<16; ++i) 
      copy_transposed_state(&current_states[i*NWORDS_STATE],
			    states,
			    i);


    
    nbytes_non_equal = memcmp(current_states, next_states, 16*HASH_STATE_SIZE);
    is_corrupt = (0 != nbytes_non_equal);

    if (is_corrupt) {
      printf("found a curropt state at step=%lu\n", step);
      fprintf(fp_verify, "found a curropt state at step=%lu\n", step);
      print_byte_txt("current", (u8*)current_states, HASH_STATE_SIZE*16);
      print_byte_txt("next", (u8*)next_states, HASH_STATE_SIZE*16);
    }
    elapsed = wtime() - elapsed;
    printf("rank=%d, step=%lu, is_corrupt=%d, elapsed=%0.2fsec\n",
	   myrank,
	   step,
	   is_corrupt,
	   elapsed);
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
  fclose(fp);




  puts("Going to check...");
  verify_middle_states(rank, size, MPI_COMM_WORLD);


  MPI_Finalize();

  return 0;
  
}
