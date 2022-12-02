// Long message attack on sha256

// this program can't attack more than 128bits
// we make the following modifications:

// TODO
// 1- store intermediate values in  sha256.c fil
// 3- check the code is sound


// define which sha256 to use 
#include "numbers_shorthands.h"
#include "sha256.h"

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

#include "config.h"
#include "timing.h"
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include <mpi.h>




void phase_i_store(const size_t n,
		   size_t server_capacity[],
		   size_t global_difficulty,
		   size_t nservers){

  // ==========================================================================+
  // Hash a long message of zeros. Store the digest in a file k where          |
  // 0<= k < nservers. To decide where to store the digest h,                  |
  //  compute k := h mod  nservers                                             |
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
  u32 ones = (1LL<<global_difficulty) - 1;

  /// timing variables
  double start = 0;
  double elapsed = 0;

  // INIT SHA256 
  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);

  // store the hash value in this variable
  
  u32 state[8] = {
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


  /* if one server gets filled, it will */
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
      fwrite(state, sizeof(u32), 3, data_to_servers[k]);
      server_capacity[k] -= 1; // We need to store less blocks now
      ++nhashes_stored;
    }
    

    // decide should not stop or should stop?
    /* not the most optimal implementation */
    should_NOT_stop = 0;
    for (size_t i=0; i<nservers; ++i) {
      /*-----------------------------------------------------------------------*/
      /* Since we reduce the server capacity by one each time when we add it to*/
      /* its file. If any server has capacity larger than zero then it means   */
      /* that we should. continue hashing till all servers capacities are met. */
      /*-----------------------------------------------------------------------*/  
      /* should_NOT_stop == 0 iff all servers_capacities are 0; */
      should_NOT_stop |= (server_capacity[i] > 0) ;
    }
    

    // + save states after required amount of intervals
    if (nhashes_stored % interval == 0) {
      FILE* states_file = fopen(states_file_name, "a");
      fwrite(state, sizeof(u32), 8, states_file);
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


  u32 digest[3] = {0};
  // Read the first two elements as one 64bit unsigned int 
  u64 store_as_idx;
  
  // Exit as soon we have exausted the file or the memory
  for (; (mem_nhashes > 0) && (fp_nhashes > 0); fp_nhashes--){
    // size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
    fread(digest, sizeof(u32), 3, fp);
    // Read the first two elements as one 64bit unsigned int 
    store_as_idx = ((u64*) digest) [0] ;
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


static inline void increment_as_128(u32 ctr[4]){
  // ===========================================================+
  // Summary: Increment ctr by one  as if it is  128bit number  |
  // -----------------------------------------------------------+
  // Probably this is an over kill                              |
  // -----------------------------------------------------------+
  
  unsigned __int128* ctr_128bit = (unsigned __int128*) ctr;
  ctr_128bit[0] += 1;  
}


//+ todo pass offset as a pointer
//+ todo use typedef unsigned __int128 u128
static void find_hash_distinguished(u32 M[16],     /* in,out*/
				    u32 Mstate[8], /* in,out*/
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
  const static u32 init_state[8] = { 0x6a09e667, 0xbb67ae85,
					  0x3c6ef372, 0xa54ff53a,
					  0x510e527f, 0x9b05688c,
					  0x1f83d9ab, 0x5be0cd19 };
  
  /* placeholder to increment 128 bits by one */
  u128* offset = (u128*) M;
  
  
  while (1) { // loop till a dist pt found
    ++(offset[0]); /* increase the first 128bits by one */
    memcpy(Mstate, init_state, 32);
    sha256_single(Mstate, (unsigned char*) M);

    /* is its digest is a distinguished pt? */
    if (Mstate[0] && dist_test == 0)
      return;
  }
}



static inline u8 to_which_server(u8 MState[NWORDS_DIGEST*WORD_SIZE],
		     u32 difficulty,
		     u32 nservers)
{
  // ==========================================================================+
  // Summary: Given a state. Decide to which server we send it to.             |   
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `Mstate` : Array of bytes of the digest.                                  |
  // `difficulty`: Number of left/right most bits that has to be zero          |
  // `nservers`: How many servers we have.                                     |
  // --------------------------------------------------------------------------+


  /* Given a state. Decide to which server we send to */
  /* """ One potential idea is to do: */
  /* server <-- h % nserver           */
  /* h'     <-- h / nserver """       */
  /* with the above suggestion it might be better to treat the digest as an */
  /* array of bytes                                                         */

  /* We assume that the difficulty is a multiple of 8. todo fix this */
  u32 difficulty_in_bytes = difficulty/8;

  /* this has to be changed is we wisht to change it to u16 */
  u8 snd_to_server =  MState[difficulty_in_bytes]; /* max 255 servers */

  return snd_to_server;
  
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




  //+ todo MPI INIT commands
  
 


  // --------------------- INIT MPI & Shared Variables ------------------------|
  int nservers, myrank;
  
  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nservers);

  /* send digests, 1 digest ≡ 256 bit */
  // int nthreads = omp_get_max_threads(); // no need for this?!




  
  if (myrank >= nservers){ /* generate hashes, send them*/
    /* I am a sending processor, I only generate hashes and send them */
    /* flatten u32 snf_buf_dgst[NSERVERS][MY_QUOTA*NWORDS_DIGEST]; */
    u32* snd_buf_dgst = (u32*) malloc(sizeof(u32*)
				      *nservers
				      *PROCESS_QUOTA
				      *NWORDS_DIGEST);

    /* flatten u32 snf_buf_offst[NSERVERS][MY_QUOTA*NWORDS_OFFSET]; */
    u32* snd_buf_offst = (u32*) malloc(sizeof(u32*)
				       *nservers
				       *PROCESS_QUOTA
				       *NWORDS_OFFSET);

    u64 nfound_potential_collisions = 0; 
    u32 M[NWORDS_INPUT]; /* random word */
    u32 Mstate[8];
    
    /* Get a random message only once */
    getrandom(M, NWORDS_INPUT*WORD_SIZE, 1);

    /* random message = random things || myrank || 64bit nonce || 64bit ctr */
    u64 ctr = 0; // = 2^64 - 1 = ff...f 
    u64 nonce;
    int server_number;
    getrandom(&nonce, sizeof(u64), 1);

    M[0] = 0; /* ctr_l */ 
    M[1] = 0;/*  ctr_h */

    /* nonce */
    M[2] = nonce;
    M[3] = (nonce>>32);

    

    u32* servers_ctr = (u32*) malloc(sizeof(u32)*nservers);
    memset(servers_ctr, 0, sizeof(u32)*nservers); /* everything = 0 */


    /* generate hashes */
    //while (needed_collisions > nfound_potential_collisions) {
    while (1) { /* yay, infinite loop */


      /* it might be better to have a nested loops */
      for (size_t i=0; i<PROCESS_QUOTA; ++i){
	
	/* Find a message that produces distinguished point */
	find_hash_distinguished(M, Mstate, difficulty_level);
	//+ decide to which server to add to? 
	server_number = to_which_server((u8*) Mstate, difficulty_level, nservers);
	snd_buf_dgst[]
	++servers_ctr[server_number];
	

	
	//+ todo check if a server snd_buf has been filled
	//+ if yes, send it immediately using buffered send 
	
	/* add it to snd_buf in location reserved for the receiver  */
	


      }
  
    }
    
    
  } else {
    /* I am a receiving processor, I only probe the dictionary */
    u32* rcv_buf = (u32*) malloc(sizeof(u32)
					   *BUFF_SIZE*NWORDS_DIGEST);
    // todo fill this 

  }
  

} // no idea where it belongs, but it silences flychek


// phase iii
//+ master server combines all files
//+ master server sort the (hash, message) according to hash
//+ master server look at the long messag file, and search for
//+ hashes in the sorted list above
