// Long message attack
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
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include "util_files.h"
#include "common.h"
// ------------------- Auxililary functions phase iii --------------------------
//+ todo complete these functions




/* 1 if  dgst1 > dgst2, -1 if dgst1<dgist2, 0 if dgst1==dgst2 */
int cmp_dgst(void const* dgst1, void const* dgst2){
  return memcmp(dgst1, dgst2, N); /* comparison order: low bytes first */
}




void load_text_file_as_u64(u64* dest, FILE* fp, size_t nlines){
  // load FILE fp data into dest, assuming fp is a text file
  // where each line as a number
  
  static const int max_len = 50; // max lenght of a line inside a file
  char tmp[max_len];
  char* endptr; // for strtoull
  size_t idx = 0;
  
  
  while (!feof(fp)) {
    memset(tmp, 0, max_len);
    fgets(tmp, max_len, fp);

    dest[idx] = strtoull(tmp, &endptr, 10);
    ++idx;
    
    if (idx >= nlines)
      return;
  }

  
}


int main(int argc, char* argv[]) /* single machine */
{
  // 1- load message candidates and hahs them
  // -- copy the hashes, and sort the copied hashes 
  // 2- load the middle states, hash the long message
  // -- after hashing a message block, query if the state is among the ordered
  //    digests, if so, use linear to find the hash index in (unordered hashes).
  //    Use this index to retrieve the message candidate.
  
  
  // ----------------------------- PART 1 ------------------------------------

  /* load messages candidates, hash them, sort them */
  FILE* fp = fopen("data/receive/messages/archive", "r");
  size_t nmsgs = (get_file_size(fp)) / HASH_INPUT_SIZE;
  


  printf("============================\n"
	 "We have %lu candidates, \n"
	 "============================\n", nmsgs);

  const WORD_TYPE state_init[NWORDS_STATE] = {HASH_INIT_STATE};
  /* We have three arrays: */
  u8* msgs  = (u8*) malloc( sizeof(u8)*nmsgs*HASH_INPUT_SIZE );
  u8* dgsts = (u8*) malloc(sizeof(u8)*nmsgs*(N) );
  u8* dgsts_orderd = (u8*) malloc(sizeof(u8)*nmsgs*(N) );

  
  fread(msgs, nmsgs, HASH_INPUT_SIZE, fp);
  fclose(fp);



  /* Hash all the candidate messages. */
  /* one thread is enough, we'll parellize it if it's a bottle neck */
  WORD_TYPE state[NWORDS_STATE];
  for (size_t i=0; i<nmsgs; ++i) {
    /* copy init state, then hash  */
    memcpy(state, state_init, NWORDS_STATE*WORD_SIZE);
    hash_single(state, &msgs[i*HASH_INPUT_SIZE]);

    assert(is_dist_state((u8*)state));

    /* get dgst in dgst */
    memcpy(&dgsts[i*N], state, N);
  }
  memcpy(state, state_init, NWORDS_STATE*WORD_SIZE);


  // ----------------------------- PART 2 ------------------------------------
  // sort msg_dgst according to the digests
  memcpy(dgsts_orderd, dgsts, nmsgs*N);
  qsort( dgsts_orderd, nmsgs, N, cmp_dgst);


  
  puts("digests unordered:");
  for (int i=0; i<5; ++i) {
    print_byte_txt("", &dgsts[i*N], N);    
  }

  puts("digests ordered:");
  for (int i=0; i<5; ++i) {
    print_byte_txt("", &dgsts_orderd[i*N], N);  
  }

  // ----------------------------- PART 3 ------------------------------------
  // Let's do it in the simple way

  
  memcpy(state, state_init, NWORDS_STATE*WORD_SIZE);

  u8 M[HASH_INPUT_SIZE] = {0};
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* increment the message by one each time */
  u32 ones = (1LL<<DIFFICULTY) - 1; 
  u8* ptr = NULL;


  size_t ctr = 0;
  print_byte_txt("state init=", (u8*)state, NWORDS_STATE*WORD_SIZE);
  for (size_t i=0; i<NHASHES;) {
        hash_single(state, M);
	msg_ctr_pt[0]++; /* Increment 64bit of M by 1 */
	ctr++;
	if((state[0] & ones) == 0){
	  i++;
	  ptr = bsearch(state, dgsts_orderd, nmsgs, N, cmp_dgst);
	  if (ptr){
	    printf("Yes at %lu\n", i);
	    
	    size_t idx = linear_search((u8*)state,
				       dgsts,
				       nmsgs,
				       N);
	    printf("at index=%lu, random msg=\n", idx);
	    print_byte_array(&msgs[idx*HASH_INPUT_SIZE], HASH_INPUT_SIZE);
	    printf("long message ctr=%llu\n", ((u64*)M)[0]);
	    print_byte_txt("hash long=", (u8*)state, N);

	    WORD_TYPE state_rnd[NWORDS_STATE] = {HASH_INIT_STATE};
	    hash_single(state_rnd, &msgs[idx*HASH_INPUT_SIZE]);
	    print_byte_txt("hash rnd =", (u8*)state_rnd, N);

	    puts("----------------------------");
	  }
	  	 
	}
  }


  

  
  free(msgs);
  free(dgsts);
  free(dgsts_orderd);


} // quit the function



