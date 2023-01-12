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



int main(int argc, char* argv[]){

  /* print_attack_information(); */

  /* // pretend we are doing phase i, (just take few hashes) */
  /* WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE}; */
  /* u8 M[HASH_INPUT_SIZE] = {0}; */
  

  /* CTR_TYPE msg_ctr = 0; */
  /* CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /\* increment the message by one each time *\/ */
  /* *msg_ctr_pt = msg_ctr; /\* update the message counter as the given input *\/ */
  /* // store the hash value in this variable */
  /* /\* WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE}; *\/ */

  /* /\* treat state as bytes (have you considered union?) *\/ */
  /* /\* u8* stream_pt = (u8*) state;  *\/ */


  /* u32 ones = (1LL<<DIFFICULTY) - 1; */

  /* printf("N=%dbytes, Difficulty=%d, ones=%d\n", N, DIFFICULTY, ones); */

  /* puts("------------------------------\n"); */
  /* int nmsgs = 0; */
  /* const static u32 ones_nservers = (1LL<<LOG2_NSERVERS) - 1; */
  /* int server = 0; */
  
  /* for (; nmsgs<20;){ */
  /*   hash_single(state, M); */
  /*   msg_ctr_pt[0]++; */




    
  /*   server =  ( (((WORD_TYPE*)state)[0] >> (DIFFICULTY)) & ones_nservers) */
  /*     % NSERVERS; */
    
  /*   if ( (state[0] & ones) == 0 && (server == 1) ){ /\* it is a distinguished point *\/ */
  /*     ++nmsgs; */
  /*     printf("msg_ctr=%llu, nserver=%d\n", msg_ctr_pt[0], server); */
  /*     puts("YAY"); */
  /*     print_char((u8*) state, 32); */
  /*     puts("--------------------------------------------\n"); */
  /*   } */
  /*  } */


  /* // we will pretend that we are 2nd server i.e. server #1  */
  
  /* puts("\n=====================================\n\n"); */
  /* puts("pretends that we are phase i"); */

  /* nmsgs = PROCESS_QUOTA; */
  /* FILE* fp = fopen("data/send/digests/1", "r"); */

  /* int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES) */
  /*                   + sizeof(CTR_TYPE); /\* |dgst| + |ctr| *\/ */

  /* u8* hashes = (u8*) malloc(one_pair_size*nmsgs); */
  /* // msg||dgst */
  /* for (int i = 0; i<nmsgs; i++){ */
  /*   printf("%d\n", i); */
  /*   fread(&hashes[one_pair_size*i + sizeof(CTR_TYPE)], */
  /* 	  sizeof(u8), */
  /* 	  N-DEFINED_BYTES, */
  /* 	  fp); */
    
  /* } */

  /* dict* d = dict_new(NSLOTS_MY_NODE); */
  /* dict_add_element_to(d, &hashes[ sizeof(CTR_TYPE)]); */

  /* dict_add_element_to(d, &hashes[7*one_pair_size + sizeof(CTR_TYPE)]); */

  /* puts("\n\n"); */
  /* int found = 0; */

  /* for (int i = 0; i<8; ++i) { */
  /*   found = dict_has_elm(d, &hashes[i*one_pair_size + sizeof(CTR_TYPE)]); */
  /*   printf("%dith was found=%d\n",i, found); */
  /*   print_char(&hashes[i*one_pair_size + sizeof(CTR_TYPE)], */
  /* 	       N-DEFINED_BYTES); */
  /*   puts("++++++++++++++++++++++++++++++++++++++++++++++++++++"); */


  /* } */

  /* free(hashes); */
  /* dict_free(d); */
  /* free(d); */

/// -------------------------------------------------------------------

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=--=--=-=-=-=--=-=-=-=--=


  /// ----------------------- INIT ---------------------------///
  /// 1- INIT numerical and bytes variables:
  /* phase will rehash the long message again in parallel */
  /* here we define how many parllel processors in phase iii */
  //size_t ncores = 14; // deadweight

  int should_NOT_stop = 1;
  /* size_t nhashes_stored = 0; */

  u32 ones = (1LL<<DIFFICULTY) - 1;



  /* Actually we have discussed that we can fix the number of nhashes per */
  /* Eventhough they have different capacities, since we allow server to  */
  /* discard excessive hasing */

  /* record the whole state after each each interval has passed */
  static const size_t interval = INTERVAL; 

  /// timing variables


  // INIT SHA256

    // -INIT: The number of Hashes each server will get
  CTR_TYPE msg_ctr = 0;
  size_t nhashes_stored = 0;
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};

  u8 M[HASH_INPUT_SIZE] = {0};
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* increment the message by one each time */

  *msg_ctr_pt = msg_ctr; /* update the message counter as the given input */
  // store the hash value in this variable
  /* WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE}; */

  /* treat state as bytes (have you considered union?) */
  u8* stream_pt = (u8*) state; 

  /// INIT FILES: server files that will be send later and state

  /* char states_file_name[FILE_NAME_MAX_LENGTH]; */
  
  printf("interval=%ld, nhashes_stored≈%ld\n", interval, nhashes_stored);





  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hash an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  printf("Going to generate hashed with DIFFICULTY=%d\n", DIFFICULTY);

  
  /* if one server gets filled, it will */
 
  while (should_NOT_stop) {
    // hash and extract n bits of the digest
    hash_single(state, M);
    msg_ctr_pt[0]++; /* Increment 64bit of M by 1 */
    /* print_char(M, 64); */

    
    
    if ( (state[0] & ones) == 0){ /* it is a distinguished point */
      print_char((u8*) state, N);
      /* Decide which server is responsible for storing this digest */

      ++nhashes_stored;

      
      // decide should not stop or should stop?
      /* not the most optimal implementation */
      should_NOT_stop = (nhashes_stored < NHASHES);
    

      // + save states after required amount of intervals

    }
    
  }






  
  

}
