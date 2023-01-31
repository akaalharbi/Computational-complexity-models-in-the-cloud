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

  print_attack_information();

  printf("Estimated memory per dictionary=%lu Bytes ≈ 2^%0.2f Bytes\n",
	 dict_memory(NSLOTS_MY_NODE),
	 log2(dict_memory(NSLOTS_MY_NODE)));


  if (argc > 1){
    float l = atof(argv[1]);
    u64 given_nslots = exp2(l);
    printf("%llu\n", given_nslots);
    
    printf("For l=%0.3f\n"
	   "Estimated memory per dictionary =%lu Bytes ≈ 2^%0.2f Bytes\n",
	   l,
	   dict_memory(given_nslots),
	   log2(dict_memory(given_nslots)));

    
  }
  printf("sizeof(dict)=%lu bytes\n", sizeof(dict));

  return 0;




       
}

