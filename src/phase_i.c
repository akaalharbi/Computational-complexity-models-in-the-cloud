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
//#include <sys/random.h> // probably deadweight getrandom(void *buffer, size_t length, 1 

// -----------------------------------------------------------------------------




// @todo rename file, and truncate 
void was_state_written_on_disk(CTR_TYPE* msg_ctr, /* ou t*/
			       size_t* nhashes_stored, /* out */
			       WORD_TYPE state[NWORDS_STATE] /* out */
			       )
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
    fp = fopen("data/counters", "w");
    fclose(fp); /* remove the counter file */
    return;
  }


  does_file_exits = access("data/counters", F_OK);
  if (does_file_exits == -1){
    // delete the states file
    puts("No counters file has been found. Starting from the beginning");
    fp = fopen("data/states", "w");
    fclose(fp); /* remove the counter file */
    return;
  }


  
  fp = fopen("data/states", "r");
  size_t nstates = get_file_size(fp)/HASH_STATE_SIZE;
  printf("Found %lu states saved\n"
	 "We will truncate the file to multiple states size if necessary\n",
	  nstates);


  

  // in case of interruption, a partial state might be recorded
  // remove the partial state. 
  truncate("data/states", nstates*HASH_STATE_SIZE);    
  truncate("data/counters", nstates*sizeof(CTR_TYPE));

  if (nstates==0){
    return;
  }

  
  /* go to the last state */
  fseek(fp,
	(nstates-1)*HASH_STATE_SIZE,
	SEEK_SET); 
  
  fread(state, WORD_SIZE, NWORDS_STATE, fp);
  puts("done with states reading");
  fclose(fp);

  
  fp = fopen("data/counters", "r");
  printf("while found %lu counters\n",
	 get_file_size(fp)/(sizeof(CTR_TYPE)));
  
  fseek(fp,
	(nstates-1)*sizeof(CTR_TYPE),
	SEEK_SET); 

  puts("opened the counter file");
  // max lenght of a line inside a file
  fread(msg_ctr, sizeof(CTR_TYPE), 1, fp);

  *nhashes_stored = nstates * INTERVAL; /* approximately how many hashes stored*/

  printf("Loaded from disk: counter=%llu\n", *msg_ctr);
}




void phase_i_store(CTR_TYPE msg_ctr,
		   size_t nhashes_stored,
		   WORD_TYPE state[NWORDS_STATE]){


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
  size_t k =  0; // server index
  int should_NOT_stop = 1;
  /* size_t nhashes_stored = 0; */

  u32 ones = (1LL<<DIFFICULTY) - 1;



  /* Actually we have discussed that we can fix the number of nhashes per */
  /* Eventhough they have different capacities, since we allow server to  */
  /* discard excessive hasing */

  /* record the whole state after each each interval has passed */
  static const size_t interval = INTERVAL; 

  /// timing variables
  double start = 0;
  double end = 0;
  double elapsed = wtime();

  // INIT SHA256 
  u8 M[HASH_INPUT_SIZE] = {0};
  /* increment the message by one each time */
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M;
  /* update the message counter as the given input */
  *msg_ctr_pt = msg_ctr;



  /* state is given as input,  treat state as bytes  */
  u8* state_u8 = (u8*) state; 

  /// INIT FILES: server files that will be send later and state
  char file_name[FILE_NAME_MAX_LENGTH]; // more thnan enough to store file name
  /* char states_file_name[FILE_NAME_MAX_LENGTH]; */
  
  FILE* data_to_servers[NSERVERS];
  FILE* states_file;
  FILE* counters_file;

  // TOUCH FILES ON THE DISK
  states_file = fopen("data/states", "a");
  counters_file = fopen("data/counters", "a");

  
  for (size_t i=0; i<NSERVERS; ++i) {
    //edit file name according to server i
    snprintf(file_name, sizeof(file_name), "data/digests/%lu", i);
    printf("file_name=%s\n", file_name);

    /* append to it if it exists, todo */
    data_to_servers[i] = fopen(file_name, "a"); 
  }


  printf("interval=%ld, nhashes_stored≈%ld\n", interval, nhashes_stored);


  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hashes an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  printf("Going to generate hashed with DIFFICULTY=%d\n", DIFFICULTY);


  start = wtime(); /* get the time  */


  /* Record the initial message and the initial ctr */
  /* Record the whole state */
  fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);
  /* Record the counter  */
  fwrite(msg_ctr_pt, sizeof(CTR_TYPE), 1, counters_file);
  /* We would like to flush the data disk as soon we have them */
  fflush(states_file);
  fflush(counters_file);

  
  /* if one server gets filled, it will */
  while (nhashes_stored < NHASHES) {
    // hash and extract n bits of the digest
    hash_single(state, M);
    msg_ctr_pt[0]++; /* Increment 64bit of M by 1 */
    
    
    if ( (state[0] & ones) == 0){ /* is it a distinguished point? */
      /* Decide which server is responsible for storing this digest */
      k = to_which_server( state_u8 );
      // = ( (state[0]>>DIFFICULTY) & ones_nservers) % NSERVERS;

      /* Recall that: */
      /* h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest   */
      fwrite(&state_u8[DEFINED_BYTES], /* start  from "the rest" see above */
	     sizeof(u8), /* smallest moving unit */
	     (N-DEFINED_BYTES), /* len( (dist_pt)|| h1 ) = DEFINED_BITS */
	     data_to_servers[k]);

      ++nhashes_stored;

      // decide should not stop or should stop?
      /* not the most optimal implementation */
      /* should_NOT_stop = (nhashes_stored < NHASHES); */
    
      // + save states after required amount of intervals
      if (nhashes_stored % interval == 0) {

	end = wtime(); /*  for progress report*/
	/* FILE* states_file = fopen(states_file_name, "a"); */
	
	/* Record the whole state */
	fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);

	/* Record the counter  */
	fwrite(msg_ctr_pt, sizeof(CTR_TYPE), 1, counters_file);


        /* We would like to flush the data disk as soon we have them */
	fflush(states_file);
	fflush(counters_file);
	for (int i=0; i<NSERVERS; ++i) 
	  fflush(data_to_servers[i]);

	
	printf("%2.4f%%, ETA %0.4fsec, write %0.2f MB/S, "
	       "#hashes≈%lu,  msg_ctr=%llu\n",
	       100 * ((float) nhashes_stored) /  NHASHES,
	       (end-start) * (NHASHES-nhashes_stored)/((float) interval),
	       interval * N / ((end - start) * 1000000.0),
	       nhashes_stored,
	       msg_ctr_pt[0]);

	start = wtime();
	
      }
    }
  }

  /* Record the last  message and the last ctr */
  /* Record the whole state */
  fwrite(state, sizeof(WORD_TYPE), NWORDS_STATE, states_file);
  /* Record the counter  */
  fwrite(msg_ctr_pt, sizeof(CTR_TYPE), 1, counters_file);
  /* We would like to flush the data disk as soon we have them */
  fflush(states_file);
  fflush(counters_file);

  
  fclose(states_file);
  fclose(counters_file);
  for (int i=0; i<NSERVERS; ++i)
    fclose(data_to_servers[i]);


  elapsed = wtime() - elapsed;
  printf("done in %fsec\n", elapsed);
}

int main(){
  /// Generate distinguished points to be stored in dictionaries
  /// inside each server.

  print_attack_information(); /* */
  puts("\n========================================\n");
  printf("Going to store %0.2f kB\n", NHASHES*N / 1000.0);
  // -INIT: The number of Hashes each server will get
  CTR_TYPE msg_ctr = 0;
  size_t nhashes_stored = 0;
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};

  /* update state and msg_ctr  if was some date befoer in the disk */
  was_state_written_on_disk(&msg_ctr, &nhashes_stored, state);

  /* continue hashing  */
  phase_i_store(msg_ctr, nhashes_stored, state);
  

}
