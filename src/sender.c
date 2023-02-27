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
#include <sys/mman.h>
#include "c_sha256_avx.h"




void find_hash_distinguished(u8 Mavx[16][HASH_INPUT_SIZE], /* in*/
			     u8  digests[16 * N]/* out*/,
			     CTR_TYPE counters[16], /* out */
			     CTR_TYPE* ctr /* in , out */,
			     int* n_dist_points /* out */)

{
  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found. Mutates M  |
  //          and resests Mstate in the process. Suitable for phase ii use only|
  // --------------------------------------------------------------------------+
  // We start with message M, computed hash(M) if it is not a distinguished    |
  // point then we change M slightly. Increment the first 64 bits of M by 1    |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - Mavx[lane][HASH_INPUT_SIZE]: place holder for messages to be hashed     |
  //      using avx. they are copies of M, however counter part  increases by  |
  //      one for each lane.                                                   |
  // - M[HASH_INPUT_SIZE] : initial value of the random message, when we found |
  //   a disitingusihed point in on of Mavx, we coput Mavx[lane] to M          |
  // - Mstate: This should be the initial state of sha256.                     |
  // - dist_test: (2^nzeros - 1), e.g. a point is distinguished if it has 3    |
  //              at the end, then the test is 2^3 - 1 = 7                     |
  // --------------------------------------------------------------------------+
  // WARNING: this function can't deal with more than 32zeros as dist_test     |
  // NOTE: we are only working on the first word                               |
  // --------------------------------------------------------------------------+
  // TODO: use hash_multiple instead                                           |
  // --------------------------------------------------------------------------+

  /* no need to construct init state with each call of the function */ 
  static u32* states_avx_ptr;

  /* copy the counter part of M */


  /* This forces us to use avx512, we should not!  */
  const REG_TYPE zero = SIMD_SETZERO_SI();
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(MASK);
  static REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  static REG_TYPE cmp_vect;
  
  static u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */
  

  
  while (1) { /* loop till a dist pt found */

    // why not use avx instruction?
    for (int i=0; i<(AVX_SIZE/WORD_SIZE_BITS); ++i) {
      ((CTR_TYPE*)Mavx[i])[0] = (*ctr)+(i+1); /* increase counter part in M by 1 */
    }

    *ctr += (AVX_SIZE/WORD_SIZE_BITS); /* update counter locally */
    

    #ifdef  __AVX512F__
    /* HASH 16 MESSAGES AT ONCE */
    states_avx_ptr = sha256_multiple_x16(Mavx);
    #endif

    #ifndef  __AVX512F__
    #ifdef    __AVX2__
    /* HASH 16 MESSAGES AT ONCE */
    states_avx_ptr = sha256_multiple_oct(Mavx);
    #endif
    #endif


    // test for distinguishedn point
    /* load the last words of digests  */
    digests_last_word = SIMD_LOAD_SI(&states_avx_ptr[(N_NWORDS_CEIL - 1) * HASH_STATE_SIZE]);
    /* is it a distinguish point? */
    
    cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);
    


    
    /* This is a bit annoying */
    #ifdef __AVX512F__ 
    cmp_mask = SIMD_CMP_EPI32(cmp_vect, zero);
    #endif

    #ifndef __AVX512F__
    cmp_mask  =_mm256_movemask_epi8( SIMD_CMP_EPI32(cmp_vect, zero) ) ;
    #endif


    if (cmp_mask) { /* found at least a distinguished point? */
      
      *n_dist_points = __builtin_popcount(cmp_mask);
      int idx = 0;
      int trailing_zeros = 0;

      for (int i=0; i<*n_dist_points; ++i){
	/* Basically get the index of the set bit in cmp_mask */
	trailing_zeros = __builtin_ctz(cmp_mask); 
	idx += trailing_zeros;
	cmp_mask = (cmp_mask >> trailing_zeros) ^ 1;
	/* update counter the ith counter */
	counters[i] = Mavx[idx][0];
	copy_transposed_digest(&digests[i*N], states_avx_ptr, idx);
      }
      return; /* we're done */
    } /* end if (cmp_mask) */
  } /* end while(1) */
}


static void regenerate_long_message_digests(u8 Mavx[16][HASH_INPUT_SIZE],
					    u8* restrict work_buf,
					    u8* restrict snd_buf,
					    int nsenders,
					    MPI_Comm inter_comm)
{
  // ==========================================================================+
  // Summary: Regenrate hashes using the states file, read specific part of the|
  //          file that only this process will read. Send digests to severs.   |
  // --------------------------------------------------------------------------+
  /* M = 64bit ctr || 64bit nonce || random value */
  /* 16 messages for avx, each message differs om the other in the counter part */
  /* except counter, they are all the same */
  u32* statesAVX;
  
  

}


static void generate_random_digests(u8 Mavx[16][HASH_INPUT_SIZE],
				    u8* restrict work_buf, /* servers digests here */
				    u8* restrict snd_buf, /* snd_buf */
				    MPI_Comm inter_comm)
{
  u32* statesAVX;
}



void sender( MPI_Comm local_comm, MPI_Comm inter_comm)
{

  // ------------------------------------------------------------------------+
  // I am a sending processor. I only generate hashes and send them.         |
  // Process Numbers: [NSERVERS + 1,  nproc]                                 |
  //-------------------------------------------------------------------------+

  // ----------------------------- PART 0 --------------------------------- //
  int myrank, nsenders;
  MPI_Comm_rank(local_comm, &myrank);
  MPI_Comm_size(local_comm, &nsenders);

  /* random message base */
  u8 M[HASH_INPUT_SIZE]; 
  u8 Mavx[16][HASH_INPUT_SIZE]; /* non-transposed! */
  u8 state_avx[16][HASH_STATE_SIZE]; /* non-transposed */



  /* |dgst| + |ctr| */
  size_t const one_pair_size = sizeof(u8)*N + sizeof(CTR_TYPE); 
  size_t const  nbytes_per_server = one_pair_size * PROCESS_QUOTA;
  /* { (server0 paris..) | (server1 pairs...) | ... | (serverK pairs...) } */
  u8* work_buf = (u8*) malloc( nbytes_per_server * NSERVERS );
  u8* snd_buf = (u8*) malloc(nbytes_per_server);

  /* /\* Decide where to place the kth digest in server i buffer *\/ */
  /* int server_number = -1; /\* what is this for ? *\/ */
  /* size_t offset = 0; /\* which index within a  server buffer should we pick *\/ */

  /* /\* pos i: How many messages we've generated to be sent to server i? *\/ */
  /* u64 servers_ctr[NSERVERS] = {0}; */
  /* int snd_ctr = 0; /\* how many messages have been sent *\/ */


  
  /* Get a random message only once */
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* counter pointer */
  getrandom(M, HASH_INPUT_SIZE, 1);
  msg_ctr_pt[0] = 0; /* zeroing the first 64bits of M */

  /* // copy the the random message to all avx messages */
  /* for (int i = 0; i<16; ++i)  */
  /*   memcpy(Mavx[i], M, HASH_INPUT_SIZE); */

  // print the template. this is not necessary.
  char txt[50];
  snprintf(txt, sizeof(txt), "sender #%d template=", myrank);
  print_byte_txt(txt, M,HASH_INPUT_SIZE);
  puts("\n");

  
  // ----------------------------- PART 1 --------------------------------- //
  // Regenrate the long message in parallel!                                //
  regenerate_long_message_digests(Mavx,
				  work_buf,
				  snd_buf,
				  nsenders,
				  inter_comm);

  // ----------------------------- PART 2 --------------------------------- //
  // 1-  Sen the initial input to all receiving servers 
  MPI_Allgather(M, HASH_INPUT_SIZE, MPI_UNSIGNED_CHAR, NULL, 0, NULL, inter_comm);


  // ------------------------------ PART 3 ----------------------------------- +
  // 1- generate disitingusihed hashes                                         |
  // 2- decide which server should probe the disitingusihed hash               |
  // 3- if server i buffer has enough hashes, send them immediately.           |
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}                               |
  // dgst := N bytes of the digest                                             |
  // ctr  := (this is a shortcut that allows us not to send the whole message) |
  //---------------------------------------------------------------------------+



  
  // ------------------------------ PART 4 ----------------------------------- +
  // Generate hashes and send them
  generate_random_digests(Mavx, work_buf, snd_buf, inter_comm);
  
  /* double time_start = wtime();  */
  /* printf("sender #%d: done init mpi, and sharing the its template." */
  /* 	 "going to generate hashes \n", myrank); */
  
  /* // find_hash_distinguished_init();  */


  /* while(1) { /\* when do we break? never! *\/ */
  /*   /\* Find a message that produces distinguished point *\/ */
  /*   find_hash_distinguished(Mavx, */
  /* 			    M, /\* save the message that generates dist here *\/ */
  /* 			    Mstate); /\* save the distinguished state here *\/ */

  /*   //+ decide to which server to add to?  */
  /*   server_number = to_which_server((u8*) Mstate); */


  /*   /\* 1st term: go to server booked memory, 2nd: location of 1st free place*\/ */
  /*   offset = server_number * nbytes_per_server */
  /*          + servers_ctr[server_number] * one_pair_size; */
  /*   // recall que one_pair_size =  |dgst| + |ctr| - |known bits| */




  /*   // record a pair (msg, dgst), msg is just the counter in our case */
  /*   /\* record the counter  *\/ */
  /*   memcpy(&snd_buf[ offset ], */
  /* 	   M, */
  /* 	   sizeof(CTR_TYPE) ); */

  /*   /\* After the counter save N-DEFINED_BYTES of MState *\/ */
  /*   memcpy( &snd_buf[offset + sizeof(CTR_TYPE)], /\* copy digest to snd_buf[offset] *\/ */
  /* 	    ((u8*)Mstate) + DEFINED_BYTES, /\* skip defined bytes, @todo skip the left most bytes *\/ */
  /* 	    N-DEFINED_BYTES ); /\* nbytes to be sent, compressed state. *\/ */

  /*   servers_ctr[server_number] += 1; */
    
    
  /*   if (servers_ctr[server_number] >= PROCESS_QUOTA){ */
  /*     printf("===============================================\n" */
  /*            "sender #%d -> recv #%d before sending %0.4fsec\n" */
  /* 	     "===============================================\n\n", */
  /* 	     myrank, server_number, wtime() -  time_start ); */

  /*     time_start = wtime(); */
  /*     MPI_Send(&snd_buf[server_number*nbytes_per_server], */
  /* 		PROCESS_QUOTA*one_pair_size, */
  /* 		MPI_UNSIGNED_CHAR, */
  /* 		server_number, */
  /* 		TAG_SND_DGST, */
  /* 		MPI_COMM_WORLD); */

  /*     /\* if (snd_ctr > 0) { *\/ */
  /*     /\* 	printf("to %d\n", server_number); *\/ */
  /*     /\* 	print_char(&snd_buf[server_number*nbytes_per_server], PROCESS_QUOTA*one_pair_size); *\/ */
  /*     /\* 	print_attack_information(); *\/ */
  /*     /\* } *\/ */

  /*     /\* else { *\/ */
  /*     /\* 	printf("i am %d to %d\n", myrank, server_number); *\/ */
  /*     /\* 	print_char(&snd_buf[server_number*nbytes_per_server], PROCESS_QUOTA*one_pair_size); *\/ */
  /*     /\* } *\/ */


  /*     ++snd_ctr; */

  /*     printf("-----------------------------------------------\n" */
  /*            "sender #%d -> recv #%d sending done %0.4fsec\n" */
  /* 	     "-----------------------------------------------\n\n", */
  /* 	     myrank, server_number, wtime() -  time_start ); */

  /*     servers_ctr[server_number] = 0; */
  /*     time_start = wtime(); */
  /*   } */
  /* } */ 


  free(snd_buf);
  free(work_buf);
  
  return; // au revoir.

} // MPI_Finalize


