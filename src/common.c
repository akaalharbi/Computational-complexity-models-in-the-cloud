// Long message attack
#include "numbers_shorthands.h"
#include "hash.h"

#include "dict.h"

// deadweight
// #include <bits/types/struct_timeval.h>

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h> // access functoin 
//#include "memory.h" // memory monitor 
//#include <sys/time.h> // moved timing.h 
//#include <assert.h>
#include "config.h"
#include "timing.h"
//#include "types.h" // probably deadweight
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include "util_files.h"
#include "c_sha256_avx.h"


void print_attack_information(){
  printf("\nL=%d, L_IN_BYTES=%d, N=%d, NHASHES=%llu,\n"
         "DIFFICULTY=%d, |idx| = %dbytes, NSEERVERS=%d,\n"
         "NSLOTS_MY_NODE=%llu, NPROBES_MAX=%d, VAL_SIZE=%d\n"
         "NDEFINED BYTES=%d, NCND_NEEDED=%llu≈2^%0.4f,\n"
	 "NDISCARDED_BITS=%d\n"
	 "AVX_SIZE=%dbits, mask_ones=%d\n"
	 "NPROBES LOOK UP = %d\n",
	 L,
	 L_IN_BYTES,
	 N,
	 NHASHES,
	 DIFFICULTY,
	 MIN(L_IN_BYTES, N-DEFINED_BYTES-VAL_SIZE_BYTES),
	 NSERVERS,
	 NSLOTS_MY_NODE,
	 NPROBES_MAX,
	 VAL_SIZE_BYTES,
	 DEFINED_BYTES,
	 NNEEDED_CND,
	 log2(NNEEDED_CND),
	 DISCARDED_BITS,
	 AVX_SIZE,
	 (1<<DIFFICULTY) - 1,
	 (int) (NPROBES_MAX/SIMD_LEN));

  puts("-------------------------------\n");
  
}



u32 to_which_server(u8 MState[HASH_STATE_SIZE])
{
  // ==========================================================================+
  // Summary: Given a state. Decide to which server we send it to.             |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `Mstate` : Array of bytes of the digest.                                  |
  // --------------------------------------------------------------------------+
  // Given a state. Decide to which server we send to                          |
  // """ One potential idea is to do:                                          |
  // server <-- h % nserver                                                    |
  // h'     <-- h / nserver """                                                |
  //
  // --------------------------------------------------------------------------+
  
  const static u32 ones_nservers = (1LL<<LOG2_NSERVERS) - 1;

  /* 1- convert MState to a WORD_TYPE (dist_bits || nserver || rest)32bits */
  /* 2- remove the distinguished bits by shifting  (nserver || rest ) */
  /* 3- keep only the bits that holds nserver (nserver) it may have extra bit */
  /* 4- Compute server number by taking computing mod nservers */
  u32 snd_to_server  = ((WORD_TYPE*)MState)[0] ;

  snd_to_server = (( snd_to_server >> DIFFICULTY) & ones_nservers) % NSERVERS;

  return snd_to_server;
}



void find_hash_distinguished_old(u8 M[HASH_INPUT_SIZE], /* in, out*/
				    WORD_TYPE Mstate[NWORDS_STATE], /* out*/
				    CTR_TYPE* ctr, /* in, out */
				    const size_t dist_test /* in */)
{
  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found. Mutates M  |
  //          and resests Mstate in the process. Suitable for phase ii use only|
  // --------------------------------------------------------------------------+
  // We start with message M, computed hash(M) if it is not a distinguished    |
  // point then we change M slightly. Increment the first 64 bits of M by 1    |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - M[16] : initial value of the random message, we hash this in the first  |
  //           in first trial, then change it to M xor ctr                     |
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
  const static WORD_TYPE init_state[NWORDS_STATE] = {HASH_INIT_STATE};
  


  /* increments the first sizeof(CTR_TYPE)*8 bits of M by 1 */

  
  while (1) { /* loop till a dist pt found */
    ++(*ctr); /* increase counter part in M by 1 */
   
    memcpy(Mstate, init_state, 32);
    /* todo  use hash multiple */
    /* figure out the number of words from config.h */
    hash_single(Mstate,  M);
    
    /* is its digest is a distinguished pt? */
    /* see 1st assumption in config.h  */
    if ( (((WORD_TYPE*) Mstate)[0] & dist_test) == 0){

      return; /* we got our distinguished digest */
    }

      
  }
}



void cpy_transposed_state(u32* Mstate, u32* tr_state,  int lane){
  for (int i = 0; i<NWORDS_STATE; ++i) {
    Mstate[i] = tr_state[lane + i*16];
  }
}


void find_hash_distinguished(u8 Mavx[16][HASH_INPUT_SIZE], /* in*/
			     u8 M[HASH_INPUT_SIZE], /* in, out */
			     WORD_TYPE Mstate[NWORDS_STATE], /* out*/
			     const size_t dist_test /* in */)
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
  CTR_TYPE local_ctr = ((CTR_TYPE*) M)[0]; 


  
  while (1) { /* loop till a dist pt found */
    
    for (int i=0; i<(AVX_SIZE/WORD_SIZE_BITS); ++i) {
      ((CTR_TYPE*)Mavx[i])[0] = local_ctr+(i+1); /* increase counter part in M by 1 */
    }

    local_ctr += (AVX_SIZE/WORD_SIZE_BITS); /* update counter locally */
    

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

    
    for (int i=0; i< (AVX_SIZE/WORD_SIZE_BITS); ++i) { // 8 for avx2,  16 for avx512
      
      if ((states_avx_ptr[i] & dist_test) == 0) {
	/* copy the message that produces distinguish point along with its counter */
	cpy_transposed_state(Mstate, states_avx_ptr, i);
	memcpy(M, Mavx[i], HASH_INPUT_SIZE);
	return;
      }
    }
    
  }
}



int is_dist_state(WORD_TYPE state[NWORDS_STATE*WORD_SIZE]){
  static const WORD_TYPE ones = (1LL<<DIFFICULTY) - 1;

  return ( (state[0] & ones) == 0 );
}


int is_dist_msg(u8 M[HASH_INPUT_SIZE]){
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};
  hash_single(state, M);

  return is_dist_state(state);
  
}
