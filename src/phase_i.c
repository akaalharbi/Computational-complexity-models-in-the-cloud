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




void was_state_written_on_disk(CTR_TYPE* msg_ctr,
			       size_t* nhashes_stored,
			       WORD_TYPE state[NWORDS_STATE])
{
  // ==========================================================================+
  // Summary: If phase_i_store was interrupted during its excution, load saved |
  //          state and counter from file. Set msg_ctr := found_counter, and   |
  //          memcpy(state, found_state, NWORDS_STATE*WORD_SIZE);              |
  // --------------------------------------------------------------------------+
  int does_file_exits = access("data/states", F_OK);
  if (does_file_exits == -1){
    puts("No states file has been found. Starting from the beginning");
    return;
  }


  does_file_exits = access("data/counters", F_OK);
  if (does_file_exits == -1){
    puts("No counters file has been found. Starting from the beginning");
    return;
  }
    

  FILE* fp = fopen("data/states", "r");
  size_t nstates = get_file_size(fp)/(NWORDS_STATE*WORD_SIZE);
  printf("Found %lu states saved\n", nstates);
  
  if (nstates==0){
    puts("");
    return;
  }
    
  
  fseek(fp,
	(nstates-1)*NWORDS_STATE*WORD_SIZE,
	SEEK_SET); 
  
  fread(state, WORD_SIZE, NWORDS_STATE, fp);
  puts("done with states reading");
  fclose(fp);

  fp = fopen("data/counters", "r");
  puts("opened the counter file");
  static const int max_len = 21; // max lenght of a line inside a file
  
  char tmp[max_len];
  char* endptr; // for strtoull
  size_t ctr_old = 0; /* Last read will be always zero, see man fgets */
  
  while (!feof(fp)) {
    *msg_ctr = ctr_old; /* update the message till we get the last ctr */
    
    memset(tmp, 0, max_len);
    fgets(tmp, max_len, fp);

    ctr_old = strtoull(tmp, &endptr, 10);
  }
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
  // 0<= k < nservers. It will store N bytes of digest (N defined in config.h) |
  // To decide which server gets the digest h, compute k := h1 mod  nservers   |
  // where h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest.     |
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
  // TODO:                                                                     |
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
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* increment the message by one each time */
  *msg_ctr_pt = msg_ctr; /* update the message counter as the given input */
  // store the hash value in this variable
  /* WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE}; */

  /* treat state as bytes (have you considered union?) */
  u8* stream_pt = (u8*) state; 

  /// INIT FILES: server files that will be send later and state
  char file_name[FILE_NAME_MAX_LENGTH]; // more thnan enough to store file name
  /* char states_file_name[FILE_NAME_MAX_LENGTH]; */
  
  FILE* data_to_servers[NSERVERS];
  FILE* states_file;
  FILE* counters_file;

  // TOUCH FILES ON THE DISK
 
  // fopen("data/states", "w"); @todo do we need this?
  /* create a string that will become a file name  */
  /* snprintf(states_file_name, */
  /* 	   sizeof(states_file_name), */
  /* 	   "data/states"); */

  states_file = fopen("data/states", "a");
  counters_file = fopen("data/counters", "a");
  /* fclose(states_file); // we will open this file again in few occasions */
  
  for (size_t i=0; i<NSERVERS; ++i) {
    //edit file name according to server i
    snprintf(file_name, sizeof(file_name), "data/send/digests/%lu", i);
    printf("file_name=%s\n", file_name);

    /* append to it if it exists, todo */
    data_to_servers[i] = fopen(file_name, "a"); 
  }

  // Init coutners before the beginning of the attack
  //interval = nhashes_stored / ncores;
  printf("interval=%ld, nhashes_stored≈%ld\n", interval, nhashes_stored);





  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hash an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  printf("Going to generate hashed with DIFFICULTY=%d\n", DIFFICULTY);
  start = wtime(); /* get the time  */
  
  /* if one server gets filled, it will */
 
  while (should_NOT_stop) {
    // hash and extract n bits of the digest
    hash_single(state, M);
    msg_ctr_pt[0]++; /* Increment 64bit of M by 1 */
    /* print_char(M, 64); */

    
    
    if ( (state[0] & ones) == 0){ /* it is a distinguished point */

      /* Decide which server is responsible for storing this digest */
      k = to_which_server((u8*) state);
	//( (state[0]>>DIFFICULTY) & ones_nservers) % NSERVERS;

      /* Recall that: */
      /* h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest   */
      fwrite(stream_pt+DEFINED_BYTES, /* start  from "the rest" see above */
	     sizeof(u8), /* smallest moving unit */
	     N-DEFINED_BYTES, /* len( (dist_pt)|| h1 ) = DEFINED_BITS */
	     data_to_servers[k]);

    

      ++nhashes_stored;


      
      // decide should not stop or should stop?
      /* not the most optimal implementation */
      should_NOT_stop = (nhashes_stored < NHASHES);
    

      // + save states after required amount of intervals
      if (nhashes_stored % interval == 0) {

	end = wtime(); /*  for progress report*/
	/* FILE* states_file = fopen(states_file_name, "a"); */
	
	/* Record the whole state */
	fwrite((WORD_TYPE*) stream_pt, sizeof(WORD_TYPE), NWORDS_STATE, states_file);
	
	/* Record the counter  */
	fprintf(counters_file, "%llu\n", msg_ctr_pt[0]);
	// We would like to flush the data disk as soon we have them
	fflush(states_file);
	fflush(counters_file);
	for (int i=0; i<NSERVERS; ++i) 
	  fflush(data_to_servers[i]);

	
	printf("%2.4f%% Generating %lu distinguished points took"
	       "%0.2fsec, nhashes_stored≈%lu\n",
	       100 * ((float) nhashes_stored) /  NHASHES,
	       interval, end - start, nhashes_stored);

	// elapsed += end - start; 
	start = wtime();
	
      }
    }
    
  }
  fclose(states_file);
  fclose(counters_file);
  for (int i=0; i<NSERVERS; ++i)
    fclose(data_to_servers[i]);


  // 
  elapsed = wtime() - elapsed;
  printf("done in %fsec\n", elapsed);
}

int main(){
  /// Generate distinguished points to be stored in dictionaries
  /// inside each server.
   
  // -INIT: The number of Hashes each server will get
  CTR_TYPE msg_ctr = 0;
  size_t nhashes_stored = 0;
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};
  was_state_written_on_disk(&msg_ctr, &nhashes_stored, state);

  // If we have not done the computations before:
  // start from the beginning
  phase_i_store(msg_ctr, nhashes_stored, state);
  print_attack_information();
}
