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

// ----------------------------------------------------------------------------


// phase iii
//+ master server has all the potential collisions messages and their digest
//+ a python script will combine them into a single file
//+ hash all message candidates
//+ master server sort the (hash, message) according to hash
//+ master server look at the long messag file, and search for
//+ hashes in the sorted list above


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
  const WORD_TYPE state_init[NWORDS_STATE] = {HASH_INIT_STATE};
  printf("============================\n"
	 "We have %lu candidates, \n"
	 "============================\n", nmsgs);
  
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

    /* get dgst in dgst */
    memcpy(&dgsts[i*N], state, N);
  }


  // ----------------------------- PART 2 ------------------------------------
  // sort msg_dgst according to the digests
  memcpy(dgsts_orderd, dgsts, nmsgs*N);
  qsort( dgsts_orderd, nmsgs, N, cmp_dgst);

  puts("digests unordered:");
  for (int i=0; i<nmsgs; ++i) {
    print_byte_txt("", &dgsts[i*N], N);    
  }

  puts("digests ordered:");
  for (int i=0; i<nmsgs; ++i) {
    print_byte_txt("", &dgsts_orderd[i*N], N);  
  }

  // ----------------------------- PART 3 ------------------------------------
  // Load middle states and counters
  
  fp = fopen("data/states", "r");

  size_t nmiddle_states  = get_file_size(fp)/(NWORDS_STATE*WORD_SIZE);

  printf("nmiddle_states=%lu, INTERVAL=%llu\n", nmiddle_states, INTERVAL);
  /* assert( nmiddle_states == INTERVAL); */

  WORD_TYPE* middle_states = (WORD_TYPE*) malloc(sizeof(WORD_TYPE)
						 *NWORDS_STATE
						 *WORD_SIZE
						 *INTERVAL);

  CTR_TYPE* counters = (CTR_TYPE*) malloc(sizeof(CTR_TYPE)*INTERVAL);

  /* Load all states: we've 1024 states, we don't getting all of them in RAM */
  fread(middle_states, INTERVAL*NWORDS_STATE, WORD_SIZE, fp);
  fclose(fp);

  /* load counters: they are written as a text file  */
  fp = fopen("data/counters", "r");
  load_text_file_as_u64(counters, fp, nmiddle_states);
  fclose(fp);
  
  /* How many state does a thread handle */
  size_t thread_load = nmiddle_states/omp_get_max_threads();
  /* size_t thread_load = nmiddle_states/1; */

  printf("----------------------------------------------------------\n"
	 "We have %lu middle states, load/thread %lu, nthreads=%d\n"
	 "-----------------------------------------------------------\n",
	 nmiddle_states, thread_load, omp_get_max_threads());
  

  double start_time = wtime();
  #pragma omp parallel
  {
    void* ptr; // result of binary search
    
    int thread_num = omp_get_thread_num();
    size_t start = thread_load*thread_num;
    size_t end = thread_load*(1+thread_num);
    u64 ctr_priv = -1;
    
    WORD_TYPE state_priv[NWORDS_STATE];
    WORD_TYPE test = (1LL<<DIFFICULTY) - 1;
    u8 msg_priv[HASH_INPUT_SIZE] = {0};
    double start_time_priv;


    if (thread_num == omp_get_max_threads() - 1)
      end = nmiddle_states; /* last thread gets chunk + nmsg % nthreads */

    printf("thd%d, [start=%lu, end=%lu)\n", thread_num, start, end);
    
    
    for (size_t i = start; i<end; ++i) {
      start_time_priv = wtime();
      
      ctr_priv = counters[i];

      // get the middle state 
      memcpy(state_priv,
	     &middle_states[i*NWORDS_STATE*WORD_SIZE],
	     NWORDS_STATE*WORD_SIZE);
      
      // between each middle states there are INTERVAL distinguished points
      // we only hash distinguished points because the h(random message) is
      // also a distinguished point.
      for (size_t j=0; j<INTERVAL; ++j) {
	find_hash_distinguished(msg_priv, state_priv, &ctr_priv, test);

	//ptr = bsearch(state_priv, dgsts_orderd, N, nmsgs, cmp_dgst);
	ptr = linear_search_ptr((u8*) state_priv, dgsts_orderd, nmsgs, N);

	if (ptr != NULL) {/* the binary search was successful! */

	  /* remeber digest if the firs N bytes of state  */
          #pragma omp critical
	  {
	    printf("found a collision at %lu\n", i);
	    printf("hash long message at %lu:\n", i);

	    size_t msg_idx = linear_search((u8*) state_priv, dgsts, nmsgs, N);
	    u8* ptr_msg_collide = &msgs[msg_idx];
	

	
	    print_byte_array((u8*) state_priv, N);
	
	    printf("while the the following message:\n");
	    print_byte_array(ptr_msg_collide, HASH_INPUT_SIZE);
	    puts("produce the following hash:");
	    print_byte_array(&dgsts[msg_idx], N);

	    
	  }// end critical region
	  
	} // end if condition


      }
      printf("thd%d done the task %lu in %0.2fsec)\n",
	     omp_get_thread_num(), i, wtime() - start_time_priv);
    } // end for loop, thread's main work
    
  } // end parallel region

  printf("it took %0.2f\n", wtime() - start_time );


  free(msgs);
  free(dgsts);
  free(dgsts_orderd);

} // quit the function



