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
#include <dirent.h>
#include <unistd.h>
// ------------------- Auxililary functions phase iii --------------------------
//+ todo complete these functions

static void truncate_messages(){
  /// Truncate all files in data/messages to a multiple of HASH_INPUT_SIZE bytes
  char file_name[FILE_NAME_MAX_LENGTH];
  size_t file_size; 
  size_t nmsgs;
  FILE* fp;

  for (size_t i = 0; i<NSERVERS; ++i){
    snprintf(file_name, FILE_NAME_MAX_LENGTH, "data/messages/%lu", i);
    fp = fopen(file_name, "r");
    file_size = get_file_size(fp);
    nmsgs = file_size/HASH_INPUT_SIZE;
    printf("file %s had %lu bytes\n", file_name, file_size);
    
    truncate(file_name, HASH_INPUT_SIZE*nmsgs);
    file_size = get_file_size(fp);
    printf("file %s has %lu bytes\n", file_name, file_size);

    fclose(fp);
  }
}




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
  /* Truncate messages to multiple of HASH_INPUT_SIZE */
  truncate_messages();
  
  /* load messages candidates, hash them, sort them */
  FILE* fp = fopen("data/messages/archive", "r");
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
  /* one thread is enough, it will be parallelized when it's a bottleneck  */
  WORD_TYPE state[NWORDS_STATE];
  for (size_t i=0; i<nmsgs; ++i) {
    /* copy init state, then hash  */
    memcpy(state, state_init, HASH_STATE_SIZE);
    hash_single(state, &msgs[i*HASH_INPUT_SIZE]);

    /* a sanity check  */
    assert(is_dist_state(state));

    /* get dgst in dgst */
    memcpy(&dgsts[i*N], state, N);
  }
  memcpy(state, state_init, HASH_STATE_SIZE);


  // ----------------------------- PART 2 ------------------------------------
  // sort msg_dgst according to the digests
  memcpy(dgsts_orderd, dgsts, nmsgs*N);
  qsort( dgsts_orderd, nmsgs, N, cmp_dgst);


  puts("digests unordered: (first 5)");
  for (int i=0; i<5; ++i) {
    print_byte_txt("", &dgsts[i*N], N);    
  }

  puts("digests ordered: (first 5)");
  for (int i=0; i<5; ++i) {
    print_byte_txt("", &dgsts_orderd[i*N], N);  
  }


  // ----------------------------- PART 3 ------------------------------------
  // First get all states: 
  FILE* fstates = fopen("data/states", "r");
  FILE* fctrs   = fopen("data/counters", "r");
  
  u64 nmiddle_states = get_file_size(fstates)/HASH_STATE_SIZE;
  /* between each middle state there are: */
  size_t nhashes_in_interval = NHASHES / nmiddle_states;
  printf("There %llu middle states\n", nmiddle_states);

  WORD_TYPE* middle_states = (WORD_TYPE*) malloc(nmiddle_states*HASH_STATE_SIZE);
  CTR_TYPE* middle_ctr = (CTR_TYPE*) malloc(nmiddle_states*sizeof(CTR_TYPE));

  fread(middle_states, WORD_SIZE, nmiddle_states*NWORDS_STATE, fstates);
  fread(middle_ctr,   sizeof(CTR_TYPE), nmiddle_states, fctrs);  

  fclose(fstates);
  fclose(fctrs);
  

  print_byte_txt("init state      ", (u8*) state_init, HASH_STATE_SIZE);
  print_byte_txt("1st middle state", (u8*) middle_states, HASH_STATE_SIZE);
  printf("1st counter = %llu\n", middle_ctr[0]);
  printf("nhashes in interval = %lu, nhashes=%llu, mod=%llu\n",
	 nhashes_in_interval, NHASHES, NHASHES % nmiddle_states);



  // ----------------------------- PART 4 ------------------------------------




  #pragma omp parallel for
  for(size_t ith_state=0; ith_state<(nmiddle_states-1); ++ith_state){
    double start = wtime();
    u8 M_priv[HASH_INPUT_SIZE] = {0};
    CTR_TYPE* M_ctr_pt_priv = (CTR_TYPE*) M_priv;
    CTR_TYPE next_ctr = middle_ctr[ith_state+1];
    WORD_TYPE state_priv[HASH_STATE_SIZE];
    u8* srearch_ptr_priv = NULL;

    

    

    memcpy(state_priv,
	   &middle_states[ith_state*NWORDS_STATE],
	   HASH_STATE_SIZE);

    M_ctr_pt_priv[0] = middle_ctr[ith_state];

    printf("ctr=%llu, ", ((u64*) M_priv)[0]);

    for (; M_ctr_pt_priv[0]<next_ctr; ) {
      hash_single(state_priv, M_priv);
      ++M_ctr_pt_priv[0];



      
      if(is_dist_state(state_priv)){
	srearch_ptr_priv = bsearch(state_priv, dgsts_orderd, nmsgs, N, cmp_dgst);

	if (srearch_ptr_priv){
	  printf("Yes at %llu\n", M_ctr_pt_priv[0]);
	    
	  size_t idx = linear_search((u8*)state_priv,
				     dgsts,
				     nmsgs,
				     N);

          printf("at index=%lu, random msg=\n", idx);
	  print_byte_array(&msgs[idx*HASH_INPUT_SIZE], HASH_INPUT_SIZE);
	  printf("long message ctr=%llu\n", ((u64*)M_priv)[0]);
	  print_byte_txt("hash long=", (u8*)state_priv, N);

	  WORD_TYPE state_rnd[NWORDS_STATE] = {HASH_INIT_STATE};
	  hash_single(state_rnd, &msgs[idx*HASH_INPUT_SIZE]);
	  print_byte_txt("hash rnd =", (u8*)state_rnd, N);

	  puts("----------------------------");
	}
	  	 
      }
      
    }

        printf("ith_state=%lu, next counter=%llu, done in %0.2fsec, thrd%d\n",
	       ith_state, next_ctr, wtime() - start, omp_get_thread_num());

    
  }
  


  free(msgs);
  free(dgsts);
  free(dgsts_orderd);
  free(middle_states);
  free(middle_ctr);
  
} // quit the function



