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
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"



// mask for distinguished point
#define DIST_MASK 0x7 // test if the last three bits are zeros

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///


// was there a cycle in PHASE I
int is_there_duplicate = 0;
int idx_cycle = -1;



void phase_i_store(const size_t n,
		   size_t server_capacity[],
		   size_t server_difficulty[],
		   size_t nservers){

  // ==================================================================================+
  // Hash a long message of zeros. Store the digest in a file k where 0<= k < nservers |
  // To decide where to store the digest h, compute k := h mod  nservers               |
  // We allow each server to adjust its own difficulty level (not sure we need that )  |
  // or its distinguished point format.                                                |
  // Compute h -> decides which server k -> check server difficulty -> decide to store |
  // h or discard it.                                                                  |
  // ----------------------------------------------------------------------------------|
  // INPUTS:                                                                           |
  // `n`: hash digest length, our goal is 96-bit                                       |
  // `server_capacity[]` : array of size nservers,  entry i contains how many blocks   |
  //                      n server i will store in its dictionary                       |
  // `server_difficulty[]` : entry i contains d s.t. server i only accepts digests h   |
  //                         h <= d                                                    |
  // `nservers` : how many servers we should prepare for                               |
  // ==================================================================================+
  




  /// ----------------------- INIT ---------------------------///
  /// 1- INIT numerical and bytes variables:
  size_t ncores = 14; // @tidy @todo get the value directly from config.h
  size_t k =  0; // server index
  int should_NOT_stop = 1;
  size_t nhashes_stored = 0; // 
  size_t interval = 1;


  /// timing variables
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  // INIT SHA256 
  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);

  // store the hash value in this variable
  uint64_t digest[2] = {0, 0};
  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };
  /// INIT FILES: server files that will be send later and state
  char file_name[40]; // more than enough to store file name
  char states_file_name[40];
  
  FILE* data_to_servers[nservers];
  FILE* states_file;
 
  // TOUCH FILES ON THE DISK
 
  fopen("data/states", "w");
  gettimeofday(&begin, 0);
  snprintf(states_file_name, sizeof(states_file_name), "data/%lu_state", begin.tv_sec);
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
  gettimeofday(&begin, 0);
  
  while (should_NOT_stop) {
    // hash and extract n bits of the digest
    sha256_single(state, M);
    truncate_state32bit_get_digest(digest, state, n);

    // Decide which server is responsible for storing this digest
    k = digest[0] % nservers;
    if (digest[1] >= server_difficulty[k]){
      /* if (k==0) { */
      /* 	printf("k=%ld, digest1,0=%016lx%016lx\n", k,digest[1], digest[0]); */
      /* } */

      // This is a distinguished point, we store maximally 128bits
      // in the distant future we may regret this decision.
      fwrite(digest, sizeof(uint64_t), 2, data_to_servers[k]);
      server_capacity[k] -= 1; // We need to store less blocks now
      ++nhashes_stored;
    }
    

    // decide should_stop or not?
    for (size_t i=0; i<nservers; ++i) {
      if (server_capacity[i] > 0){
	should_NOT_stop = 1;
	break; // from the inner loop
      }
    }
    // @todo start from here
    // + save states after required amount of intervals
    if (nhashes_stored % interval == 0) {
      FILE* states_file = fopen(states_file_name, "a");
      fwrite(state, sizeof(uint32_t), 8, states_file);
      // We would like to flush the data disk as soon we have them
      fclose(states_file);
    }
    
  }

  /// TIMING record PHASE I time elapsed ///
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  ///
  /// write it in a file  ///
  printf("done in %fsec, ", elapsed);
  /// -------------------///

}



void long_message_attack(const size_t n_of_bits, const double l, FILE* fp){


  

  
  /*  Mount long message attack on truncated sha256 */
  /// INPUT:
  ///  - n_of_bits: compression functions output size in bits
  ///  - l: 2^l = how many message blocks we need to store their intermediate values
  ///  -fp: a pointer to file where we will store the benchmark information (time,)
  ///      it will write in the first line it encounters then 
  /// todo nblocks instead of n_of_blocks!

  /// PROCESS THE PRARAMETERS
  is_there_duplicate = 0; // global variable to detect cycles
  size_t n_of_blocks = (size_t) ceil(pow(2.0, l));
  // We will use this to reconstruct the long message again
  size_t ncores = 14; // @tidy @todo get the value directly from config.h 
  
  
  // each thread will create number of registers @todo adapte the formula below
  double memory_estimate = (64  + 32)*omp_get_max_threads() + dict_memory(n_of_blocks); // edit me


  /// write it in a file fp
  fprintf(fp, "%lu, %d, NAN, %0.2fkb, ",
	  n_of_bits, (int) l, memory_estimate);


  printf("Memory estimates %0.2fkb, \n", memory_estimate);
  /// ----------------- IF VERBOSE ENABLED ------------------ ///
  #ifdef VERBOSE_LEVEL
  printf("Memory estimates %0.2fkb, \n", memory_estimate);
  #endif // VERBOSE_LEVEL
  /// --------------- END IF VERBOSE ENABLED ---------------- ///



  
  /// ----------------------- INIT ---------------------------///
  /// timing
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;
  gettimeofday(&begin, 0);

  
  // store values in a dictionary
  // dict_new(#elemnets will be stored, element's size)
  // assuming all element have the same size
  dict* d = dict_new(n_of_blocks); 

  // store the states of the long message after some defined interval
  FILE* fstate = fopen("data/states", "wb");
  fclose(fstate); // for now we only need to clean the file 
  /// -------------- PHASE I ------------------------  ///
  /// First phase hash an extremely long message
  // M0 M1 ... M_{2^l}
  // Mi entry will evaluated on the fly


  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  uint64_t digest[2] = {0, 0};
  // INIT SHA256 

  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };


  // hash a long message (for now it's a series of zeros)
  // number of messages recorded 
  size_t nmsg_rec = 0; // e.g. how many distinguished points we have seen so far 
  // when nmsg_rec mod interval == 0; record the state in a file
  size_t interval = n_of_blocks / ncores;
  
  while (nmsg_rec < n_of_blocks){
    // number of (distinguished messages == n_of_blocks)
    sha256_single(state, M);

    // get the digest from the state
    truncate_state32bit_get_digest(digest, state, n_of_bits);
    
    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///

    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

    if (nmsg_rec % interval == 0){
      //+ @todo
      //+ record (state) at somefile
      // should we close and open the file each time the if condition
      // becomes true? I think yes, we get the panelty on a local laptop 14 times only.
      FILE* fstate = fopen("data/states", "a");
      fwrite(state, sizeof(uint32_t), 8, fstate);
      fclose(fstate);
    }
      
    // @todo do we need to record the value in a dictionary? Or should we do it later?
    dict_add_element_to(d, digest);
    nmsg_rec++;

    /// ----------------- IF VERBOSE ENABLED ------------------ ///
    #ifdef VERBOSE_LEVEL
    printf("value=%lu\n", nmsg_rec);
    printf("digest=0x%016lx%016lx\n", digest[1], digest[0]);
    puts("state=");
    print_char((char*)state, 32);
    puts("\n");
    #endif // VERBOSE_LEVEL
    /// --------------- END IF VERBOSE ENABLED ---------------- ///
  } ///             DONE WITH HASHING THE LONG MESSAG            ///

  /// ----------------- IF VERBOSE ENABLED ------------------ ///

  #ifdef VERBOSE_LEVEL
  puts("dictionary after hash");
  dict_print(d);
  puts("-----------------");
  #endif // VERBOSE_LEVEL
  /// --------------- END IF VERBOSE ENABLED ---------------- ///


  /// TIMING record PHASE I time elapsed ///
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  ///
  /// write it in a file  ///
  fprintf(fp, "%fsec, ", elapsed);
  /// -------------------///

  /// record memeory usage ///
  long res_mem  = 0;
  long virt_mem = 0;
  get_memory_usage_kb(&res_mem, &virt_mem);
  fprintf(fp, "RES=%lukb, VmSize=%lukb, ", res_mem, virt_mem);
  ///------------------------------------------------------///


  /// ------------------ PHASE II -----------------------  ///
  // Second phase: hash a random message till a collision is found
  /// check if we have a duplicate first
  
  int collision_found = 0; // shared by all threads and dict.c file
  size_t ctr = 0; // how many random messages we've tried
  size_t idx = 0; // the index of the message to be removed
  BYTE random_message[64] = {0};

  /// The code now never checks for duplicates :(
  /* if (is_there_duplicate){ */
  /*   // this corresponds to having cycles while hashing the long message */
  /*   // i.e. h(mi  ... mj) = h(mi ... mk) where k>i */
  /*   fprintf(fp, "%f, 0, %d cycle\n", elapsed, idx_cycle); */
  /*   collision_found = 1; */
  /*   dict_free(d); */
  /*   free(d); */
  /*   return; */
  /* } */
  

  /// ------------------- ///

 
  // parallel search which are independent of each other
  #ifdef _OPENMP
  printf("Phase II: we will use %d threads, n=%lu, l=%d\n",
	 omp_get_max_threads(), n_of_bits, (int) l);
  // printf("Max nthreads=%d\n", omp_get_max_threads());

  
  #endif // _OPENMP

  /// ---------  PARELLLEL SEARCH ----------- ///
  #pragma omp parallel shared(collision_found)
  {
    //printf("Thread %02d\n", omp_get_thread_num());


    uint32_t state_init_priv[8] = {
      0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };
    //+ todo add more structured search 
    BYTE random_message_priv[64] = {0};
    uint64_t *random_message_priv64 = (uint64_t*) random_message_priv;
    random_message_priv64[0] = omp_get_thread_num()*( (1LL<<63)/((uint64_t)omp_get_max_threads()) );
    /* printf("the interval will b divided as %d intervals, each interval os size %llu\n", omp_get_max_threads(),  (1LL<<63)/(omp_get_max_threads()) ); */
    /* printf("thread %02d with offset %lu, as normal arith=%llu\n", omp_get_thread_num(), random_message_priv64[0], omp_get_thread_num()*( (1LL<<63)/(omp_get_max_threads()) )); */
    
    uint64_t digest_priv[2] = {0};
    uint32_t state_priv[8] = {0};
    size_t idx_priv = 0;
    size_t ctr_priv = 0;



    
    while (!collision_found) {
      // get a new different message
      random_message_priv64[0] += 1ULL;

      // clean previously used values
      memcpy(state_priv, state_init_priv, sizeof(state_init_priv));

      // hash
      sha256_single(state_priv, random_message_priv);

      // extract the results
      truncate_state32bit_get_digest(digest_priv, state_priv, n_of_bits);

      /// ------------ DISTINGUISHED POINTS (if enabled) ------------- ///
       /// If distinguished points feature was enabled  during compile ///
      /// time. 
      #ifdef DISTINGUISHED_POINTS
      // we skip hashes
      if ( (digest_priv[0]&DIST_MASK) != 0) 
	continue; // skip this element
      #endif

      //+ @todo this section needs more writing
      // test for collision and print the results if successful.
      // idx_priv := 0 if digest_priv doesn't exist in the dictionary
      idx_priv = dict_get_value(d, digest_priv);
      //if (dict_has_key(d, digest_priv) ){ // maybe extra condition && !collision_found is needed
      if (idx_priv ){ //

	int does_it_collide = 1; //collides_at(random_message_priv, n_of_bits, idx_priv);
	if (does_it_collide){
          #pragma omp critical
	  { // update a shared values
	    collision_found = 1;
	    idx = idx_priv; 
	    memcpy(random_message, random_message_priv, 64);

            #ifdef VERBOSE_LEVEL
	    printf("digest_rand=0x%016lx%016lx\n", digest_priv[1], digest_priv[0]);
	    printf("ctr_priv=%lu\n", ctr_priv);
	    print_char((char*)random_message_priv, 64);
	    puts("---end---\n\n");
           #endif
	  }
	  // collision found
	  break; // we don't care about the rest of the loop
	}
      }
      ++ctr_priv;
    }

    // record the number of trials
    #pragma omp critical
    {
      ctr += ctr_priv;
    }
  }


  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  /// record PHASE I time elapsed ///
  /// write it in a file  ///
  fprintf(fp, "%fsec, %lu, %lu\n", elapsed, ctr, idx); // last writing

  /// -------------------///

  // record the message
  char file_name[40];
  snprintf(file_name, sizeof(file_name), "data/messages/%lu_%d",  n_of_bits,  (int) l);
  FILE* fm = fopen(file_name, "w");
  fwrite(random_message, 1, 64, fm);
  fclose(fm);

  // Free memeory
  dict_free(d);
  //free(d);


}




int main(int argc, char* argv []){
  // attack(size_t n_of_blocks, size_t n_of_bits)

  puts("starting with storing message");
  #define nservers 16
  size_t servers_capacity[nservers];
  size_t server_difficulty[nservers];
  for (int i=0; i<nservers; ++i) {
    servers_capacity[i] = 1LL<<25;
    server_difficulty[i] = 0;
  }
  phase_i_store(96, servers_capacity, servers_capacity, nservers);
  puts("done with registering the message");

  
  
  int n = 0;
  float l = 0;

  char first_line[] =   "n l estimated_time estimated_memory time_phase_i real_memory_usage virt_memory_usage time_all_attack trials idx_message\n";

  char directory_name[] = "data/stats/";




  /// How to handle the arguments given by the command line:
  /// case 1: no arguments, for oarsub we will assume now it will
  ///         be same as case 5 with nmax:=100, nmin:=71, lmin:=32, lmax:=32
  /// case 2: n is given /// I think we should remove this case 
  /// case 3: n l are given
  /// case 5: nmax nmin lmax lmin are given



  
  if (argc == 3){ // ./long_message_attack n l
  // get the input from the user
  n = atoi(argv[1]);
  l = atof(argv[2]);

  // variable file name
  char file_name[40];
  snprintf(file_name, sizeof(file_name), "%s%d_%d_stats.txt",
	   directory_name, (int) n, (int) l);
  printf("filename=%s\n", file_name);

  /// todo write the value of n and l explicitly in the file 
  //  FILE* fp = fopen("statistics/n_l_stats.txt", "w");
  FILE* fp = fopen(file_name, "w");
  fprintf(fp, "%s", first_line);
	    

  long_message_attack(n, l, fp);

  fclose(fp);

  return 0;

  }  else if (argc == 5) { // ./long_message_attack nmax nmin lmax lmin

    // supply n_max n_min l_max l_min
    int n_max = atoi(argv[1]);
    int n_min = atoi(argv[2]);
    int l_max = atoi(argv[3]);
    int l_min = atoi(argv[4]);
    // l = atof(argv[2]);

    printf("n_max=%d,  n_min=%d, l_max=%d, l_min=%d\n",
	   n_max, n_min, l_max, l_min);
    // variable file name
    char file_name[50];
    snprintf(file_name, sizeof(file_name), "%s%d_%d_%d_%d_stats.txt", directory_name, n_max, n_min, l_max, l_min);
    printf("filename=%s\n", file_name);



	     
    FILE* fp = fopen(file_name, "w");
    fprintf(fp, "%s", first_line);    fclose(fp);
    
    /// loop over n_min <= n <= n_max, l_min <= l <= l_max
    for (int n1=n_min; n1 <= n_max; ++n1){
      // when l > n/2 then the we expect to be a cylce in phase I
      for (int l=l_min; (l<=l_max && l<= (n1>>1) ); ++l){
	// opening file multiple time to results as soon we as we have it
	FILE* fp = fopen(file_name, "a");
	long_message_attack(n1, l, fp);
	fclose(fp);
	// puts("");
      }
    }


  }
  
  else {
    
    
    puts("==================================");
    puts("Welcome to the Long message attack");
    puts("USAGE: ./long_message_attack n_max n_min l_max l_min\n"
	 "This will test all paris (n, l) such that:\n"
	 "n_min<= n <= n_max, l_min <= l <= l_max\n"
	 "The results will be saved in the file:\n"
	 "nmax_nmin_lmax_lmin_stats.txt in the statistics folder.");
    puts("USAGE: ./long_message_attack n l\n"
	 "This will only record the usage of n l."
	 "The corresponding satistics will be found in"
	 "n_l_stats.txt");
    
    puts("USAGE: ./long_message_attack n\n"
	 "This will try 0<n and increase n by 1 each time."
	 "It will save the statistics in the file satistics/stats.txt");
    puts("n: 0<n<257, the number of bits in the output of the compression function\n"
	 "l:positive double, 2^l is the number of blocks");
    return 0;
  }


  
  return 0;
}



