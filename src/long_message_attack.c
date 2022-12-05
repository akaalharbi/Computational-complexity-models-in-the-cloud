// Long message attack on sha256

// this program can't attack more than 128bits
// we make the following modifications:

// TODO
// 1- store intermediate values in  sha256.c fil
// 3- check the code is sound


// define which sha256 to use 
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

#include "config.h"
#include "timing.h"
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include <mpi.h>


//---------------------------- UTILITY FUNCTIONS -------------------------------


static void find_hash_distinguished(WORD_TYPE M[NWORDS_INPUT], /* in, out*/
				    WORD_TYPE Mstate[NWORDS_STATE], /* out*/
				    const size_t dist_test /* in */)
{
  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found. Mutates M  |
  //          and resests Mstate in the process. Suitable for phase ii use only|
  // --------------------------------------------------------------------------+
  // We start with message M, computed hash(M) if it is not a distinguished    |
  // point then we change M slightly. Increment the first 64 bits of M by 1    |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - M[16] : initial value of the random message, we hash this in the first  |
  //           in first trial, then change it to M xor ctr                     |
  // - Mstate: This should be the initial state of sha256.                     |
  // - dist_test: (2^nzeros - 1), e.g. a point is distinguished if it has 3    |
  //              at the end, then the test is 2^3 - 1 = 7                     |
  // --------------------------------------------------------------------------+
  // WARNING: this function can't deal with more than 32zeros as dist_test     |
  // NOTE: we are only working on the first word                               |
  // --------------------------------------------------------------------------+
  // TODO: Don't use hash_multiple instead                                     |
  // --------------------------------------------------------------------------+

  /* no need to construct init state with each call of the function */ 
  const static WORD_TYPE init_state[NWORDS_STATE] = {HASH_INIT_STATE};
  


  /* increments the first sizeof(CTR_TYPE)*8 bits of M by 1 */
  CTR_TYPE* ctr_pt = (CTR_TYPE*) M;
  
  while (1) { /* loop till a dist pt found */
    ++(ctr_pt[0]);
   
    memcpy(Mstate, init_state, 32);
    hash_single(Mstate, (u8*) M);
    
    /* is its digest is a distinguished pt? */
    /* see 1st assumption in config.h  */
    if ( ((u32*) Mstate)[0] && dist_test == 0)
      return; /* we got our distinguished digest */
  }
}

// -----------------------------------------------------------------------------

static inline u32 to_which_server(u8 MState[NWORDS_DIGEST*WORD_SIZE])
{
  // ==========================================================================+
  // Summary: Given a state. Decide to which server we send it to.             |   
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `Mstate` : Array of bytes of the digest.                                  |
  // --------------------------------------------------------------------------+
  // Given a state. Decide to which server we send to                          |
  // """ One potential idea is to do:                                          |
  // server <-- h % nserver                                                    |
  // h'     <-- h / nserver """                                                |
  //
  // --------------------------------------------------------------------------+
  
  const static u32 ones_nservers = (1LL<<LOG2_NSERVERS) - 1;
  /* 1- convert MState to a WORD_TYPE (dist_bits || nserver || rest)32bits */
  /* 2- remove the distinguished bits by shifting  (nserver || rest ) */
  /* 3- keep only the bits that holds nserver (nserver) it may have extra bit */
  /* 4- Compute server number by taking computing mod nservers */
  u32 snd_to_server  = ( (((WORD_TYPE*)MState)[0] >> DIFFICULTY)
			 & ones_nservers) % NSERVERS;

  return snd_to_server;
}



void phase_i_store(size_t server_capacity[]){
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
  // `server_capacity[]` : array of size nservers,  entry i contains how many  |
  //                       blocks server i will store in its dictionary        |
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
  size_t ncores = 14; /* here we define how many parllel processors in phase iii */
  size_t k =  0; // server index
  int should_NOT_stop = 1;
  size_t nhashes_stored = 0; // 

  u32 ones = (1LL<<DIFFICULTY) - 1;



  /* Actually we have discussed that we can fix the number of nhashes per */
  /* Eventhough they have different capacities, since we allow server to  */
  /* discard excessive hasing */
  size_t nhashes=0;
  for (size_t i=0; i<NSERVERS; ++i) {
    /* this should coincide with number of hashes in config.h */
    nhashes += server_capacity[i]; 
  }

  /* record the whole state after each each interval has passed */
  size_t interval = nhashes>>10; 
  
  /// timing variables
  double start = 0;
  double elapsed = 0;

  // INIT SHA256 
  u8 M[NWORDS_INPUT*WORD_SIZE] = {0};

  // store the hash value in this variable
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};

  /* treat state as bytes (have you considered union?) */
  u8* stream_pt = (u8*) state; 

  /// INIT FILES: server files that will be send later and state
  char file_name[40]; // more than enough to store file name
  char states_file_name[40];
  
  FILE* data_to_servers[NSERVERS];
  FILE* states_file;
 
  // TOUCH FILES ON THE DISK
 
  fopen("data/states", "w");
  snprintf(states_file_name,
	   sizeof(states_file_name),
	   "data/%llu_state",
	   (u64) wtime());

  states_file = fopen(states_file_name, "w");
  fclose(states_file); // we will open this file again in few occasions
  
  for (size_t i=0; i<NSERVERS; ++i) {
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
    hash_single(state, M);
    

    if ((state[0] & ones) == 0){ /* it is a distinguished point */
      
      /* Decide which server is responsible for storing this digest */
      k = to_which_server((u8*) state);
	//( (state[0]>>DIFFICULTY) & ones_nservers) % NSERVERS;

      /* Recall that: */
      /* h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest   */
      fwrite(stream_pt+DEFINED_BYTES, /* start  from "the rest" see above */
	     sizeof(u8), /* smallest moving unit */
	     N-DEFINED_BYTES, /* len( (dist_pt)|| h1 ) = DEFINED_BITS */
	     data_to_servers[k]);
      
      server_capacity[k] -= 1; // We need to store less blocks now
      ++nhashes_stored;

      
      // decide should not stop or should stop?
      /* not the most optimal implementation */
      should_NOT_stop = 0;
      for (size_t i=0; i<NSERVERS; ++i) {
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
	/* Record the whole state */
	fwrite((WORD_TYPE*) stream_pt, sizeof(WORD_TYPE), NWORDS_STATE, states_file);


	// We would like to flush the data disk as soon we have them
	fclose(states_file);
      }
    }
    
  }

  /// TIMING record PHASE I time elapsed ///
  elapsed = wtime() - start;
  /// write it in a file  ///
  printf("done in %fsec, ", elapsed);
  /// -------------------///


}

//+ edit from here
void phase_i_load(dict *d, size_t fp_nhashes, FILE *fp)
{
  // ==========================================================================+
  // Summary: Load hashes from the file *fp, and store them in dict *d         |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d`: Dictionary that will keep elements from *fp.                        |
  // `*fp` : File contain number of hashes larger than nelements               |
  // `fp_nhashes`  : The nhashes available in the file.                        |
  // `mem_nhashes` : The max number of hashes we are allowed to store in this  |
  //                   machine. The dictionary *d may not accepts all elements.|
  // ==========================================================================+


  /* Check that file exists, the file comes from external resources */  
  if (!fp){
    puts("I have been given a file that lives in nowhere");
    return; // raise an error instead
  }

  u8 stream_pt[N-DEFINED_BYTES];
  /* add as many hashes as possible */
  for (;; fp_nhashes--){
    fread(stream_pt, sizeof(u8), N-DEFINED_BYTES, fp);
    /* it adds the hash iff nprobes <= NPROBES_MAX */
    dict_add_element_to(d, stream_pt); 
  }
  fclose(fp); // close the file
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
