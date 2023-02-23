// Long message attack
#include "numbers_shorthands.h"
#include "hash.h"

#include "dict.h"
#include "common.h"

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
//#include <sys/random.h> // probably deadweight getrandom(void *buffer, size_t
//length, 1

// @todo we can deduce counter in was_state_written_on_diesk from the nstates
// we should simplify the code, and having few break points.


// -----------------------------------------------------------------------------


/* static void truncate_digests(){ */
/*   /// Truncate all files in data/digests to a multiple of N bytes */
/*   char file_name[FILE_NAME_MAX_LENGTH]; */
/*   size_t file_size;  */
/*   size_t ndigests; */
/*   FILE* fp; */

/*   for (size_t i = 0; i<NSERVERS; ++i){ */
/*     snprintf(file_name, FILE_NAME_MAX_LENGTH, "data/digests/%lu", i); */
/*     fp = fopen(file_name, "r"); */
/*     file_size = get_file_size(fp); */
/*     ndigests = file_size/N; */
/*     printf("file %s had %lu bytes\n", file_name, file_size); */
    
/*     truncate(file_name, N*ndigests); */
/*     file_size = get_file_size(fp); */
/*     printf("file %s has %lu bytes\n", file_name, file_size); */

/*     fclose(fp); */
/*   } */
/* } */


// @todo rename file, and truncate 
int load_checkpoint(WORD_TYPE state[NWORDS_STATE], /* out */
		    CTR_TYPE* msg_ctr)
{
  // ==========================================================================+
  // Summary: If phase_i_store was interrupted during its excution, load saved |
  //          state and counter from file. Set msg_ctr := found_counter, and   |
  //          memcpy(state, found_state, NWORDS_STATE*WORD_SIZE);              |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - *msg_ctr: 
  FILE* fp;
  
  int does_file_exits = access("data/states", F_OK);
  if (does_file_exits == -1){
    puts("No states file has been found. Starting from the beginning");
    return 0;
  }

  
  fp = fopen("data/states", "r");
  size_t nstates = get_file_size(fp)/HASH_STATE_SIZE;
  printf("Found %lu states saved\n"
	 "We will truncate the file to multiple states size if necessary\n",
	  nstates);

  if (nstates==0)
    return 0;


  
  /* go to the last state */
  fseek(fp,
	(nstates-1)*HASH_STATE_SIZE,
	SEEK_SET); 
  
  fread(state, WORD_SIZE, NWORDS_STATE, fp);
  
  puts("state loaded = ");
  print_char((u8*) state, HASH_STATE_SIZE);
  fclose(fp);

  // Ensure that digests are multiple of N
  truncate("data/states", nstates*HASH_STATE_SIZE);

  *msg_ctr = (nstates - 1) * INTERVAL; /* approximately how many hashes stored*/


  /* truncate_digests(); */
  printf("Loaded from disk: counter=%llu\n", *msg_ctr);

  return 1;
}




void phase_i_store(CTR_TYPE msg_ctr,
		   u64 nhashes_computed,
		   WORD_TYPE state[NWORDS_STATE],
		   int was_there_data){


  		   /* const size_t n, */
		   /* size_t global_difficulty, */
		   /* size_t nservers */
                   /* parameters above moved to config.h */
  // ==========================================================================+
  // Summary:                                                                  |
  // Hash a long message of zeros. Store the digest, h, in a file k where      |
  // 0<= k < NSERVERS. It will store (N-#defined_bytes) bytes of digest.       |
  // To decide which server gets the digest h, compute k := h1 mod  nservers   |
  // where h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest.     |
  // dist_pt means the bits are picked as distinguished point test             |
  // --------------------------------------------------------------------------+
  // INPUTS: CAPITAL LETTERS input are defined in config.h                     |
  //                                                                           |
  // `msg_ctr` : We hash a long message where each block content is zero except|
  //             a small part that indicates the block number.                 |
  // `state` : Pass the base state, since each block hashing change the digest | 
  // `DIFFICULTY` : Number of bits that are 0 in the first word A:=state[0]    |
  //                       i.e. is state[0]&(2**global_difficulty - 1) == 0?   |
  // `NSERVERS` : how many servers we should prepare for                       |
  // --------------------------------------------------------------------------+
  // NOTE: Bits corresponding to distinguished point and server number are not |
  //       stored in the file.                                                 |
  // NOTE: endinaness: u32 A[2] = {x, y} then (uint64*) A = { (y<<32) | x) }   |
  // --------------------------------------------------------------------------+
  // ==========================================================================+
  



  /// ----------------------- INIT ---------------------------///
  /// 1- INIT numerical and bytes variables:
  /* phase will rehash the long message again in parallel */
  /* here we define how many parllel processors in phase iii */
  //size_t ncores = 14; // deadweight
  /* size_t nhashes_stored = 0; */



  /* Actually we have discussed that we can fix the number of nhashes per */
  /* Eventhough they have different capacities, since we allow server to  */
  /* discard excessive hasing */

  /* record the whole state after each each interval has passed */
  

  /// timing variables
  double start = 0;
  double end = 0;
  double elapsed = wtime();
  u64 MB = 0;
  
  // INIT SHA256 
  u8 M[HASH_INPUT_SIZE] = {0};
  /* increment the message by one each time */
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M;
  /* update the message counter as the given input */
  *msg_ctr_pt = msg_ctr;


  /// INIT FILES: server files that will be send later and state
  FILE* states_file;


  // TOUCH FILES ON THE DISK
  states_file = fopen("data/states", "a");


  printf("interval=%lld, nhashes_stored≈%lld\n", INTERVAL, nhashes_computed);
  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hashes an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  printf("Going to generate hashes with length %u bits\n", N );


  start = wtime(); /* get the time  */

  if (!was_there_data) { /* record the first state iff there was no state file */
    /* This ensure that last state in run_i won't be duplicated in run_{i+1} */
    /* Record the whole state */
    fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);
    /* We would like to flush the data disk as soon we have them */
    fflush(states_file);
  }


  
  /* if one server gets filled, it will */
  while (1) {
    // hash and extract n bits of the digest
    hash_single(state, M);
    msg_ctr_pt[0]++; /* Increment 64bit of M by 1 */
    ++nhashes_computed;


    // + save states after required amount of intervals
    if (nhashes_computed % INTERVAL == 0) {
      printf("nhases=%llu, interval=%llu, nhashes mode interval=%llu\n",
	     nhashes_computed, INTERVAL, nhashes_computed % INTERVAL);
      end = wtime(); /*  for progress report*/
      /* FILE* states_file = fopen(states_file_name, "a"); */
	
      /* Record the whole state */
      fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);


      /* We would like to flush the data disk as soon we have them */
      fflush(states_file);

      MB = (INTERVAL * N) / (u64) 1000000;
      printf("2^%2.4f hashes, elapsed %0.2fsec,  write %0.2f MB/S\n"
	     " #hashes≈%llu, 2^%0.3f hashes/sec,  msg_ctr=%llu\n"
	     "---------------------------------\n",
	     log2(nhashes_computed),
	     (end-start),
	     MB / ((end - start) ),
	     nhashes_computed,
	     log2( INTERVAL / (u64) (end-start)),
	     msg_ctr_pt[0]);
      
      start = wtime();
    }
  }
  

  /* Record the last  message and the last ctr */
  /* Record the whole state */
  fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);
  /* Record the counter  */

  /* We would like to flush the data disk as soon we have them */
  fflush(states_file);

  fclose(states_file);

  elapsed = wtime() - elapsed;
  printf("done in %fsec\n", elapsed);
}

int main(){
  /// Generate distinguished points to be stored in dictionaries
  /// inside each server.

  print_attack_information(); /* */
  puts("\n========================================\n");
  /* printf("Going to store %0.2f kB\n", NHASHES*N / 1000.0); */
  // -INIT: The number of Hashes each server will get
  CTR_TYPE msg_ctr = 0;
  size_t nhashes_stored = 0;
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};

  /* update state and msg_ctr  if was some date befoer in the disk */
  int was_there_data = load_checkpoint( state, &msg_ctr );

  /* continue hashing  */
  phase_i_store(msg_ctr, nhashes_stored, state, was_there_data);
  

}
