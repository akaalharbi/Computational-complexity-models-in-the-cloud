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
//#include <sys/random.h> // probably deadweight getrandom(void *buffer, size_t length, 1



int main(int argc, char* argv[]){

  
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};
  u8 M[HASH_INPUT_SIZE] = {0};

  CTR_TYPE msg_ctr = 0;
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* increment the message by one each time */
  *msg_ctr_pt = msg_ctr; /* update the message counter as the given input */
  // store the hash value in this variable
  /* WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE}; */

  /* treat state as bytes (have you considered union?) */
  /* u8* stream_pt = (u8*) state;  */


  u32 ones = (1LL<<DIFFICULTY) - 1;

  printf("N=%dbytes, Difficulty=%d, ones=%d\n", N, DIFFICULTY, ones);

  puts("------------------------------\n");
  int nmsgs = 0;
  const static u32 ones_nservers = (1LL<<LOG2_NSERVERS) - 1;
  
  for (; nmsgs<20;){
    hash_single(state, M);
    msg_ctr_pt[0]++;




    

    if ( (state[0] & ones) == 0){ /* it is a distinguished point */
      ++nmsgs;
      printf("msg_ctr=%llu, nserver=%d\n", msg_ctr_pt[0],
	     ( (((WORD_TYPE*)state)[0] >> (DIFFICULTY)) & ones_nservers) % NSERVERS
);
      puts("YAY");
      print_char((u8*) state, 32);
      puts("--------------------------------------------\n");
    }
   }
}
