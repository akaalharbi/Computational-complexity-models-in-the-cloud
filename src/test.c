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
//#include <sys/random.h> // probably deadweight getrandom(void *buffer, size_t length, 1

/* 1 if  dgst1 > dgst2, -1 if dgst1<dgist2, 0 if dgst1==dgst2 */
int cmp_dgst(void const* dgst1, void const* dgst2){
  return memcmp(dgst1, dgst2, N-DEFINED_BYTES); /* comparison order: low bytes first */
}



int main(int argc, char* argv[]){
  u8 rnd_msg[HASH_INPUT_SIZE] = {0x41, 0x4d, 0x20, 0xd9, 0x00, 0x00, 0x00, 0x00, 0x93, 0x2f, 0x3d, 0x43, 0xa5, 0x15, 0x9b, 0x35, 0x50, 0x0b, 0x8d, 0xa4, 0x8a, 0xaf, 0x47, 0xb7, 0x0e, 0x3c, 0x0b, 0x0f, 0x32, 0x7b, 0xad, 0xd2, 0xfc, 0x2c, 0xfa, 0xa2, 0xe3, 0x1f, 0x19, 0x3b, 0x94, 0x75, 0x80, 0xdf, 0x1d, 0x0c, 0xa0, 0x92, 0x39, 0x72, 0xe3, 0xde, 0x29, 0x11, 0x21, 0x6c, 0x09, 0x21, 0x0c, 0xad, 0x39, 0x1f, 0x55, 0x50 };

  u64 ctr = 323907923;

  u8 M[HASH_INPUT_SIZE] = {0};
  u64* msg_ctr_pt = (u64*) M;
  
  WORD_TYPE state_rnd_msg[NWORDS_STATE] = {HASH_INIT_STATE};
  WORD_TYPE state_long_msg[NWORDS_STATE] = {HASH_INIT_STATE};

  hash_single(state_rnd_msg, rnd_msg);

  for (size_t i = 0; i<ctr-10; ++i) {
    msg_ctr_pt[0]++;
      hash_single(state_long_msg, M);

      
  }

  for (size_t i = 0; i<11; ++i) {
      msg_ctr_pt[0]++;
      hash_single(state_long_msg, M);

      printf("ctr=%llu\n",  msg_ctr_pt[0]++);
      print_byte_txt("h(lng_msg)", (u8*) state_long_msg, N);
  }

    
  print_byte_txt("h(rnd_msg)", (u8*) state_rnd_msg, N);

  return 0;




       
}

