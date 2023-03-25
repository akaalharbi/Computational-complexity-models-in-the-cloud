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

#include "common.h"

void print_attack_information(){
  printf("\nL=%f, L_RECEIVER=%f, N=%d, NHASHES=%llu,\n"
         "DIFFICULTY=%d, |idx| = %dbytes, NELELEMNTS_BUCKET=%d, NSEERVERS=%d,\n"
         "NSLOTS_MY_NODE=%llu, NPROBES_MAX=%d, VAL_SIZE=%d\n"
         "NDEFINED BYTES=%d, NCND_NEEDED=%llu≈2^%0.4f,\n"
	 "NDISCARDED_BITS=%d\n"
	 "AVX_SIZE=%dbits, mask_ones=%d\n"
	 "NPROBES LOOK UP = %d\n",
	 L,
	 L_RECEIVER,
	 N,
	 NHASHES,
	 DIFFICULTY,
	 (int) ceil((L_RECEIVER - log2(SIMD_LEN))/8.0),
	 SIMD_LEN,
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


int is_dist_state(u8 state[HASH_STATE_SIZE]){
  /* check if the last X bits  if we read digest as little endian are zeros  */
  /* terms on the left get all 1s, term on right move 1s to the end */
  static const u8 ones = ( (1LL<<DIFFICULTY) - 1) << (8-DIFFICULTY);
  u8 last_8bits = ( (u8*) &state[N-1] ) [0];
    
  /* For now, we are not going to have difficulyt more than 16 bits. It will  */
  /* take 330 years using my laptop! The longest recorded life is 122 years!  */
  return ( (  last_8bits  & ones) == 0 );
}

int is_dist_digest(u8 state[N]){
  /* check if the last X bits  if we read digest as little endian are zeros  */

  /* terms on the left get all 1s, term on right move 1s to the end */
  static const u8 ones = ( (1LL<<DIFFICULTY) - 1) << (8-DIFFICULTY);
  u8 last_8bits = ( (u8*) &state[N-1] ) [0];
    
  /* For now, we are not going to have difficulyt more than 16 bits. It will  */
  /* take 330 years using my laptop! The longest recorded life is 122 years!  */
  return ( (  last_8bits  & ones) == 0 );
}


int is_dist_msg(u8 M[HASH_INPUT_SIZE]){
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};
  hash_single(state, M);

  return is_dist_state((u8*) state);
}



u32 to_which_server(u8 state[HASH_STATE_SIZE])
{
  // ==========================================================================+
  // Summary: Given a state. Decide to which server we send it to.             |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                     nn              |
  // `Mstate` : Array of bytes of the digest.                                  |
  // --------------------------------------------------------------------------+
  // Given a state. Decide to which server we send to                          |
  // """ One potential idea is to do:                                          |
  // server <-- h % nserver                                                    |
  // h'     <-- h / nserver """                                                |
  //
  // --------------------------------------------------------------------------+
  /* 1- convert sttate to a u32 (digest || nserver || dist_pt)  */
  /*  -> (nserver || dist_pt) 32 bits */
  /* 2- remove the distinguished bits by shifting  (nserver ) */
  /* 3- keep only the bits that holds nserver (nserver) it may have extra bit */
  /* 4- Compute server number by taking computing mod nservers */
  u32 snd_to_server  = ((u32*) &state[N-4])[0];
			
  snd_to_server = ( snd_to_server >> DIFFICULTY)  % NSERVERS;

  return snd_to_server;
} 



void transpose_state(u32 dest[restrict 16*8],
		     u32 src[restrict 16*8])
{
  for (int lane = 0; lane < 16; lane++) {
    dest[lane + 0*16] = src[8*lane + 0];
    dest[lane + 1*16] = src[8*lane + 1];
    dest[lane + 2*16] = src[8*lane + 2];
    dest[lane + 3*16] = src[8*lane + 3];
    dest[lane + 4*16] = src[8*lane + 4];
    dest[lane + 5*16] = src[8*lane + 5];
    dest[lane + 6*16] = src[8*lane + 6];
    dest[lane + 7*16] = src[8*lane + 7];
  }

}


void untranspose_state(u32 dest[restrict 16*8],
		       u32 src[restrict  8*16])
{
  /* looks stupid but I am not passing dimension as an argument */
    for (int lane = 0; lane < 16; lane++) {
    dest[8*lane + 0] = src[lane + 0*16];
    dest[8*lane + 1] = src[lane + 1*16];
    dest[8*lane + 2] = src[lane + 2*16];
    dest[8*lane + 3] = src[lane + 3*16];
    dest[8*lane + 4] = src[lane + 4*16];
    dest[8*lane + 5] = src[lane + 5*16];
    dest[8*lane + 6] = src[lane + 6*16];
    dest[8*lane + 7] = src[lane + 7*16];
  }

}







void copy_transposed_digest(u8 *digest, u32 *tr_state, int lane)
{
  for (int i = 0; i<(N/WORD_SIZE); ++i) 
    memcpy(&digest[i*WORD_SIZE],
	   &tr_state[lane + i*16],
	   WORD_SIZE);


  /* copy the rest of the bytes */
  memcpy(&digest[(N/WORD_SIZE)*WORD_SIZE],
	 &tr_state[lane + (N/WORD_SIZE)*16],
	 N - (N/WORD_SIZE)*WORD_SIZE );

}




