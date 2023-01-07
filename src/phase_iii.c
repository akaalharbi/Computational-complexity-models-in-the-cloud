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

// ------------------- Auxililary functions phase iii --------------------------
//+ todo complete these functions

/* 1 if  dgst1 > dgst2, -1 if dgst1<dgist2, 0 if dgst1==dgst2 */
int cmp_dgst(void const* dgst1, void const* dgst2){
  return memcmp(dgst1, dgst2, N); /* comparison order: low bytes first */
}

/* return index of key if it is found, -1 otherwise*/
size_t linear_search(u8 *key, u8 *array, size_t array_len, size_t key_len)
{
  for (size_t i=0; i<array_len; i += key_len) {
    if ( 0 == memcmp(key, &array[i], key_len) )
      return i;
  }
  
  return -1; /* not found */
}


void print_byte_array(u8* array, size_t nbytes)
{
  for (size_t i=0; i<nbytes; ++i) 
    printf("0x%02x, ",  array[i]);
  puts("");
}
// ----------------------------------------------------------------------------


// phase iii
//+ master server has all the potential collisions messages and their digest
//+ a python script will combine them into a single file
//+ hash all message candidates
//+ master server sort the (hash, message) according to hash
//+ master server look at the long messag file, and search for
//+ hashes in the sorted list above





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
  size_t nmsgs = get_file_size(fp)/HASH_INPUT_SIZE;
  const WORD_TYPE state_init[NWORDS_STATE] = {HASH_INIT_STATE};
  printf("we have %lu candidates\n", nmsgs);
  
  /* We have three arrays: */
  u8* msgs  = (u8*) malloc( sizeof(u8)*nmsgs*HASH_INPUT_SIZE );
  u8* dgsts = (u8*) malloc(sizeof(u8)*nmsgs*(N) );

  u8* dgsts_orderd = (u8*) malloc(sizeof(u8)*nmsgs*(N) );
  
  fread(msgs, nmsgs, HASH_INPUT_SIZE, fp);
  fclose(fp);




  /* one thread is enough, we'll parellize it if it's a bottle neck */
  WORD_TYPE state[NWORDS_STATE];
  for (size_t i=0; i<nmsgs; ++i) {
    /* copy init state, then hash  */
    memcpy(state, state_init, NWORDS_STATE*WORD_SIZE);
    hash_single(state, &msgs[i*HASH_INPUT_SIZE]);

     /* get dgst in dgst */
    memcpy(&dgsts[i*N], state, N);
  }


  // ----------------------------- PART 3 ------------------------------------
  // sort msg_dgst according to the digests
  memcpy(dgsts_orderd, dgsts, nmsgs*N);
  qsort( dgsts_orderd, nmsgs, N, cmp_dgst);
  
  // ----------------------------- PART 2 ------------------------------------
  // hash the long message with each hashing probe 
  fp = fopen("data/states", "r");
  size_t nmiddle_states  = get_file_size(fp)/(NWORDS_STATE*WORD_SIZE);
  /* How many state does a thread handle */
  size_t thread_load = nmiddle_states/omp_get_max_threads();
  

  #pragma omp parallel
  {
    void* ptr; // result of binary search
    int thread_num = omp_get_thread_num();
    size_t start = thread_load*thread_num;
    size_t end = thread_load*(1+thread_num);
    if (thread_num == omp_get_max_threads() - 1)
      end = nmiddle_states -1; /* last thread gets chunk + nmsg % nthreads */
    //+ todo how to incorporate counter?
    
    WORD_TYPE state_priv[NWORDS_STATE] = {HASH_INIT_STATE};
    u8 msg_priv[HASH_INPUT_SIZE] = {0};
    
    for (size_t i = start; i<end; ++i) {
      memcpy(state_priv, state_init, NWORDS_STATE*WORD_SIZE);
      // use counter to increment msg_priv
      hash_single(state_priv, msg_priv);
      ptr = bsearch(state_priv, dgsts_orderd, N, nmsgs, cmp_dgst);

      if (ptr) {
        /* remeber digest if the firs N bytes of state  */
	#pragma omp critical
	{
	
	size_t msg_idx = linear_search((u8*) state_priv, dgsts, nmsgs, N);
	u8* ptr_msg_collide = &msgs[msg_idx];
	

	
        printf("found a collision at %lu\n", i);
	printf("hash long message at %lu:\n", i);
	print_byte_array((u8*) state_priv, N);
	
	printf("while the the following message:\n");
	print_byte_array(ptr_msg_collide, HASH_INPUT_SIZE);
	puts("produce the following hash:");
	print_byte_array(&dgsts[msg_idx], N);
	
	}// end critical region
      } // end if condition
      
    
    } // end for loop, thread's main work

  } // end parallel region
  free(fp);
  free(msgs);
  free(dgsts);
  free(dgsts_orderd);
} // quit the function



