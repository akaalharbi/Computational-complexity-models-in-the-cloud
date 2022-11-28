// Long message attack on sha256

// this program can't attack more than 128bits
// we make the following modifications:

// TODO
// 1- store intermediate values in  sha256.c fil
// 3- check the code is sound


// define which sha256 to use 
#include "sha256.h"

#include "dict.h"
#include <bits/types/struct_timeval.h>
// #include <endian.h> // @todo do we need it?
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <sys/time.h>
#include <omp.h>

#include "config.h"
#include "timing.h"
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)






// mask for distinguished point
#define DIST_MASK 0x7 // test if the last three bits are zeros

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///








void phase_i_store(const size_t n,
		   size_t server_capacity[],
		   size_t global_difficulty,
		   size_t nservers){

  // ==================================================================================+
  // Hash a long message of zeros. Store the digest in a file k where 0<= k < nservers |
  // To decide where to store the digest h, compute k := h mod  nservers               |
  // We allow each server to adjust its own difficulty level (not sure we need that )  |
  // or its distinguished point format.                                                |
  // Compute h -> decides which server k -> check server difficulty -> decide to store |
  // h or discard it.                                                                  |
  // ----------------------------------------------------------------------------------+
  // INPUTS:                                                                           |
  // `n`: hash digest length, our goal is 96-bit                                       |
  // `server_capacity[]` : array of size nservers,  entry i contains how many blocks   |
  //                      n server i will store in its dictionary                      |
  // `global_difficulty` : Number of bits that are 0 in the first word A:=state[0]     |
  //                       i.e. is state[0]&(2**global_difficulty - 1) == 0 or not?    |
  // `nservers` : how many servers we should prepare for                               |
  // NOTE: endinaness: u32 A[2] = {x, y} then (uint64*) A = { (x<<32) | y) }           |
  // ----------------------------------------------------------------------------------+
  // TODO:                                                                             |
  // - Load server capacities from a file.                                             |
  // ==================================================================================+
  




  /// ----------------------- INIT ---------------------------///
  /// 1- INIT numerical and bytes variables:
  size_t ncores = 14; // @tidy @todo get the value directly from config.h
  size_t k =  0; // server index
  int should_NOT_stop = 1;
  size_t nhashes_stored = 0; // 
  size_t interval = 1;
  uint32_t ones = (1LL<<global_difficulty) - 1;

  /// timing variables
  double start = 0;
  double elapsed = 0;

  // INIT SHA256 
  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);

  // store the hash value in this variable
  
  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };
  // we can read digest directly from the state

  /// INIT FILES: server files that will be send later and state
  char file_name[40]; // more than enough to store file name
  char states_file_name[40];
  
  FILE* data_to_servers[nservers];
  FILE* states_file;
 
  // TOUCH FILES ON THE DISK
 
  fopen("data/states", "w");
  snprintf(states_file_name, sizeof(states_file_name), "data/%lu_state", (size_t) wtime());
  states_file = fopen(states_file_name, "w");
  fclose(states_file); // we will open this file again in few occasions
  
  for (size_t i=0; i<nservers; ++i) {
    //edit file name according to server i
    snprintf(file_name, sizeof(file_name), "data/upload/%lu", i);
    printf("file_name=%s\n", file_name);

    data_to_servers[i] = fopen(file_name, "w");
    nhashes_stored += server_capacity[i];
  }

  // Init coutners before the beginning of the attack
  interval = nhashes_stored / ncores;
  printf("interval=%ld, nhashes_stores=%ld, ncores=%ld\n",
	 interval, nhashes_stored, ncores);
  nhashes_stored = 0; // we have not recorded any hash yet





  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hash an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  start = wtime();


  
  while (should_NOT_stop) {
    // hash and extract n bits of the digest
    sha256_single(state, M);
    

    // Decide which server is responsible for storing this digest

    k = (state[0]>>global_difficulty) % nservers;
    if ((state[0] & ones) == 0){ // it is a distinguished point 
      /* if (k==0) { */
      /* 	printf("k=%ld, digest1,0=%016lx%016lx\n", k,digest[1], digest[0]); */
      /* } */

      // This is a distinguished point, we store maximally 128bits
      // in the distant future we may regret this decision.
      fwrite(state, sizeof(uint32_t), 3, data_to_servers[k]);
      server_capacity[k] -= 1; // We need to store less blocks now
      ++nhashes_stored;
    }
    

    // decide should not stop or should stop?
    for (size_t i=0; i<nservers; ++i) {
      /*-----------------------------------------------------------------------*/
      /* Since we reduce the server capacity by one each time we add it to its */
      /* files. If any server has capacity larger than zero then it means that */
      /* we should. continue hashing till all servers capacities are met.      */
      /*-----------------------------------------------------------------------*/  
      if (server_capacity[i] > 0){
	should_NOT_stop = 1;
	break; // from the inner loop
      }
    }
    

    // + save states after required amount of intervals
    if (nhashes_stored % interval == 0) {
      FILE* states_file = fopen(states_file_name, "a");
      fwrite(state, sizeof(uint32_t), 8, states_file);
      // We would like to flush the data disk as soon we have them
      fclose(states_file);
    }
    
  }

  /// TIMING record PHASE I time elapsed ///
  elapsed = wtime() - start;
  ///
  /// write it in a file  ///
  printf("done in %fsec, ", elapsed);
  /// -------------------///

}


void phase_i_load(dict *d,
		  FILE *fp,
		  size_t fp_nhashes,
		  size_t mem_nhashes,
                  size_t difficulty) {

  // ===========================================================================+
  // Summary: Load hashes from the file *fp, and store the in dict *d           |
  // ---------------------------------------------------------------------------+
  // INPUTS:                                                                    |
  // `*d`: Dictionary that will keep elements from *fp.                         |
  // `*fp` : File contain number of hashes larger than nelements                |
  // `fp_nhashes`  : The nhashes available in the file.                         |
  // `mem_nhashes` : The max number of hashes we are allowed to store in this   |
  //                   machine. The dictionary *d may not accepts all elements. |
  // ---------------------------------------------------------------------------+
  // TODO:                                                                      |
  // ===========================================================================+

  /// Load file of hashes to dictionary d
  /// this is local function, it doesn't know which server it's
  //+ check available memory
  //+ available memeory >= nelements
  //++ add all elemnets of file to d
  //+ available memeory < nelements
  //++ Fix some distinguish points or store first elemnets till
  //   memeory is full

  // Check that file exists
  if (!fp){
    puts("I have been given a file that lives in nowhere");
    return; // raise an error instead
  }


  uint32_t digest[3] = {0};
  // Read the first two elements as one 64bit unsigned int 
  uint64_t store_as_idx;
  
  // Exit as soon we have exausted the file or the memory
  for (; (mem_nhashes > 0) && (fp_nhashes > 0); fp_nhashes--){
    // size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
    fread(digest, sizeof(uint32_t), 3, fp);
    // Read the first two elements as one 64bit unsigned int 
    store_as_idx = ((uint64_t*) digest) [0] ;
    // remove the distinguished points. They are shared between all 
    store_as_idx = store_as_idx >> difficulty;
    //+ should I also remove the server number?
    //+ I think this will help with clustering

    // Add element to the dictionary. If it is successful,
    // reduce mem_nhashes by one.
    mem_nhashes -= dict_add_element_to(d, store_as_idx, digest[2]);
    
  }
  fclose(fp); // close the file
}


static inline void increment_as_128(uint32_t ctr[4]){
  // ===========================================================+
  // Summary: Increment ctr by one  as if it is  128bit number  |
  // -----------------------------------------------------------+
  // Probably this is an over kill                              |
  // -----------------------------------------------------------+
  
  unsigned __int128* ctr_128bit = (unsigned __int128*) ctr;
  ctr_128bit[0] += 1;  
}

void find_hash_distinguished(uint32_t M[16],     /* in,out*/
                             uint32_t Mstate[8], /* in,out*/
                             const size_t dist_test /* in */)
{
  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found.            |
  //          Mutate M and Mstate in the process                               |
  // --------------------------------------------------------------------------+
  // We start with message M, computed sha256(M) if it is not a distinguished  |
  // point then we change M slightly. We keep a counter ctr, that gets         |
  // increased by 1 with each tiral. Compute sha256(M xor ctr).                |
  //                                                                           |
  // INPUTS:                                                                   |
  // - M[16] : initial value of the random message, we hash this in the first  |
  //           in first trial, then change it to M xor ctr                     |
  // - Mstate: This should be the initial state of sha256.                     |
  // - dist_test: (2^nzeros - 1), e.g. a point is distinguished if it has 3    |
  //              at the end, then the test is 2^3 - 1 = 7                     |
  // --------------------------------------------------------------------------+
  // WARNING: this function can't deal with more than 32zeros as dist_test     |
  // WARNING: we are only working on the first word                            |
  // --------------------------------------------------------------------------+

  // no need to construct init state with each call of the 
  const static uint32_t init_state[8] = { 0x6a09e667, 0xbb67ae85,
					  0x3c6ef372, 0xa54ff53a,
					  0x510e527f, 0x9b05688c,
					  0x1f83d9ab, 0x5be0cd19 };
  
  // uint64_t ctr = 0; // make this 128bits, maybe using c++
  static uint32_t ctr[4] = {0};
  
  while (1) { // loop till a dist pt found
    // hash message
    memcpy(Mstate, init_state, 32);
    sha256_single(Mstate, (unsigned char*) M);

    // is its digest is a distinguished pt?
    if (Mstate[0] && dist_test == 0)
      return;

    // we wish to increase the first 128-bit by 1
    increment_as_128(ctr);
    // Message xor ctr
    M[0] ^= ctr[0];
    M[1] ^= ctr[1];
    M[2] ^= ctr[2];
    M[3] ^= ctr[3];
  }
}




void phase_ii(dict* d,
	      FILE* fp_possible_collisions,
	      size_t needed_collisions,
	      size_t difficulty_level)
{ 
  // some extra arguments to communicate with other
  
  // ==========================================================================+
  // Summary: Hash many different messages till we found enough collisions     |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d` : dictionary that has hashes from phase i                            |
  // `*fp_possible_collisions`: store all heashes h that dictionary tells they |
  //                            might exists. Also, store the corresponding    |
  //                            message.                                       |
  // NOTE: (64bytes message || 32bytes hash) (64bytes message || 32bytes hash) |
  // `needed_collisions` : How many messages we need that return positive hits |
  //                       when probe the dictionary.                          |
  //---------------------------------------------------------------------------+
  // TODO:                                                                     |
  // - extra arguments to deal with other servers                              |
  // ==========================================================================+


  // ---------------------- -PARALLEL SEARCH ---------------------------|

   #pragma omp parallel
  {
    //                     INIT PRIVATE VARAIABLES                     //
    int in_dict_priv = 0;
    //+ call
    //- @todo do we need this? it is already on find_hash_distinguished
    uint32_t state_priv[8] = { 0x6a09e667, 0xbb67ae85,
			       0x3c6ef372, 0xa54ff53a,
			       0x510e527f, 0x9b05688c,
			       0x1f83d9ab, 0x5be0cd19 };

    
    // containers for hash and random message that have hash
    // satisfies the difficulty level
    uint32_t M_priv[16]; // 512-bits 
    getrandom(M_priv, 64, 1);
    // that have hash as
    uint32_t M_state_priv[8];
    uint64_t store_as_idx_priv; // first two words of M_state



    while (needed_collisions > 0) {
      // Find message that produces distinguished point
      find_hash_distinguished(M_priv, M_state_priv, difficulty_level);

      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
      // DECIDES HOW TO DEAL WITH THE HASH:
      // 1- ADD TO SERVER i BUFFER, i NEED TO BE COMPUTED
      // 2- PROBE THE LOCAL DICTIONARY
      // --------
      // I am not sure, if probing and sending should be seperated from this
      // function.
      

      // Imagine we will probe it locally:
      store_as_idx_priv = ((uint64_t*) M_state_priv)[0];
      
      in_dict_priv = dict_get_value(d, store_as_idx_priv, state_priv[2]);
      
      if (in_dict_priv) {
       #pragma omp critical
	{
	needed_collisions--;
	//+ write hash and message to file
	// todo start from here
      } // end if 
    } // end while (needed_collisions > 0)
  } // end parallel search
    //+ send the file to the master server
}  // end function
} // no idea where it belongs, but it silences flychek


// phase iii
//+ master server combines all files
//+ master server sort the (hash, message) according to hash
//+ master server look at the long messag file, and search for
//+ hashes in the sorted list above
