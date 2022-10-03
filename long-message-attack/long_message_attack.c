// Long message attack on sha256
// this program can't attack more than 128bits
// we make the following modifications:

// TODO
// 1- store intermediate values in  sha256.c fil
// 2- doxygene style
// 3- check the code is sound



#include "sha256.h"
#include "dict.h"
#include <bits/types/struct_timeval.h>
#include <endian.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "types.h"
#include "util/util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "util/memory.h"
#include <sys/time.h>
#include <omp.h>

// mask for distinguished point
#define DIST_MASK 0x7 // test if the last three bits are zeros

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///




// int (*functionPtr)(int,int);
// add_element_to_dictionary(dict *dictionary, char *key, size_t value, size_t input_size)



// was there a cycle in PHASE I
int is_there_duplicate = 0;
int idx_cycle = -1;

void nothing(dict *dictionary, dict_key* key, size_t value, size_t key_size){
  // literally do nothing;
}


void long_message_attack(size_t n_of_bits, double l, FILE* fp){
//(size_t n_of_bits, size_t n_of_blocks){    
  /*  Mount long message attack on truncated sha256 */
  /// INPUT:
  ///  - n: compression functions output size in bits
  ///  - l: 2^l = how many message blocks we need to store their intermediate values
  ///  -fp: a pointer to file where we will store the benchmark information (time,)
  ///      it will write in the first line it encounters then 
  /// todo nblocks instead of n_of_blocks!


  /// PROCESS THE PRARAMETERS
  is_there_duplicate = 0; // global variable to detect cycles
  size_t n_of_blocks = (size_t) ceil(pow(2.0, l));
  //  size_t nelements = (size_t) ceil(pow(2.0, n_of_bits));
  // int n_of_bytes = (int) ceil( (float)n_of_bits/8);
  
  // size of  long_message (lazy evaluation) + dict 
  double memory_estimate = 64 + sizeof(dict);
         memory_estimate += 2*n_of_blocks*sizeof(slot); // dictionary size
	 memory_estimate = memory_estimate / 1000.0; //kb

  /// write it in a filep
  fprintf(fp, "%lu, %d, NAN, %0.2fkb, ",
	  n_of_bits, (int) l, memory_estimate);

  printf("Memory estimates %0.2fkb, \n", memory_estimate);
  /// ----------------- IF VERBOSE ENABLED ------------------ ///
  #ifdef VERBOSE_LEVEL
  printf("Memory estimates %0.2fkb, \n", memory_estimate);
  #endif // VERBOSE_LEVEL
  /// --------------- END IF VERBOSE ENABLED ---------------- ///



  // we need to take into account also the various pointers in this program
  // printf("ESTIMATED MEMEORY %lukb\n", memory_estimate/1000);
  
  
  /// ------- INIT ----------------------------------///
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
  
  

  /// -------------- PHASE I ------------------------  ///
  /// First phase hash an extremely long message
  // M0 M1 ... M_{2^l}
  // Mi entry will evaluated on the fly


  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  SHA256_CTX ctx;
  sha256_init(&ctx);
  // we extract n_bits digest from ctx.state
  // digest is the output of the compression function after
  // truncation.
  // unsigned char digest[n_of_bytes];

  // hash a long message (for now it's a series of zeros)
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_transform(&ctx, M, n_of_bits);
    // get the digest from the state
    truncate_state_get_digest(digest, &ctx, n_of_bits);
    // add it to the dicitonary
    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///
    /// time. 
    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

    dict_add_element_to(d, digest, i);


    /// ----------------- IF VERBOSE ENABLED ------------------ ///
    #ifdef VERBOSE_LEVEL
    printf("value=%lu\n", i);
    printf("digest=0x%016lx%016lx\n", digest[1], digest[0]);
    puts("ctx.stat=");
    print_char((char*)ctx.state, 32);
    puts("\n");
    #endif // VERBOSE_LEVEL
    /// --------------- END IF VERBOSE ENABLED ---------------- ///

    
  }

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
  ///---------------------///


  /// -------------- PHASE II --------------------  ///
  // Second phase: hash a random message till a collision is found
  /// check if we have a duplicate first
  
  int collision_found = 0; // shared by all threads
  size_t ctr = 0; // how many random messages we've tried
  size_t idx = 0; // the index of the message to be removed
  BYTE random_message[64] = {0};

  if (is_there_duplicate){
    // this corresponds to having cycles while hashing the long message
    // i.e. h(mi  ... mj) = h(mi ... mk) where k>i
    fprintf(fp, "%f, 0, %d cycle\n", elapsed, idx_cycle);
    collision_found = 1;
    free(d->slots);
    free(d);
    return;
  }
  

  /// ------------------- ///

 
  // parallel search which are independent of each other
  // omp_set_num_threads( omp_get_max_threads() );
  #ifdef _OPENMP
  printf("Phase II: we will use %d threads, n=%lu, l=%d\n",
	 omp_get_max_threads(), n_of_bits, (int) l);
  // printf("Max nthreads=%d\n", omp_get_max_threads());
  #endif // _OPENMP

  #pragma omp parallel num_threads( omp_get_max_threads() )
  {
    // each variable inside is private to each thread
    BYTE random_message_priv[64] = {0};
    uint64_t digest_priv[2] = {0, 0};
    SHA256_CTX ctx_priv;
    size_t idx_priv = 0;
    //  STATE intermediate; // union // @remove
    // we will zero the excessive zeros, don't edit
    sha256_init(&ctx_priv);
    size_t ctr_priv = 0;


  
    while (!collision_found) {
      // create a random message of 64 bytes
      fill_radom_byte_array(random_message_priv, 64);


      // hash
      sha256_init(&ctx_priv);
      sha256_transform(&ctx_priv, random_message_priv, n_of_bits);
      // extract the results
      truncate_state_get_digest(digest_priv, &ctx_priv, n_of_bits);

      /// ------------ DISTINGUISHED POINTS (if enabled) ------------- ///
      /// If distinguished points feature was enabled  during compile ///
      /// time. 
      #ifdef DISTINGUISHED_POINTS
      // we skip hashes
      if ( (digest_priv[0]&DIST_MASK) != 0) 
	continue; // skip this element
      #endif

      
      // test for collision and print the results if successful.
      // idx_priv := 0 if digest_priv doesn't exist in the dictionary
      idx_priv = dict_get_value(d, digest_priv);
      //if (dict_has_key(d, digest_priv) ){ // maybe extra condition && !collision_found is needed
      if (idx_priv ){ // 
	#pragma omp critical
	{ // update a shared values
	  collision_found = 1;
	  idx = idx_priv; 
	  memcpy(random_message, random_message_priv, 64);
	}


	/* puts("Found a collision with the following details:"); */
	/* printf("#random message trials=%lu, index=%lu, M=",ctr, idx); */
	/* dict_key* intermediate = (dict_key *) ctx2.state; */
	/* print_char(intermediate->bytes, n_of_bytes); */
	/* puts(""); */
	break; // we don't care about the rest of the loop
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
  char file_name[15];
  snprintf(file_name, sizeof(file_name), "messages/%lu_%d",  n_of_bits,  (int) l);
  FILE* fm = fopen(file_name, "w");
  fwrite(random_message, 1, 64, fm);
  fclose(fm);

  // Free memeory
  free(d->slots);
  free(d);

}




int main(int argc, char* argv []){
  // attack(size_t n_of_blocks, size_t n_of_bits)
  
  int n = 0;
  float l = 0;

  char first_line[] =   "n l estimated_time estimated_memory time_phase_i real_memory_usage virt_memory_usage time_all_attack trials idx_message\n";

  #ifdef _OPENMP
  char directory_name[] = "statistics_parallel/";
  #else
  char directory_name[] = "statistics/";
  #endif


  /// How to handle the arguments given by the command line:
  /// case 1: no arguments, for oarsub we will assume now it will
  ///         be same as case 5 with nmax:=100, nmin:=71, lmin:=32, lmax:=32
  /// case 2: n is given /// I think we should remove this case 
  /// case 3: n l are given
  /// case 5: nmax nmin lmax lmin are given


  if (argc == 1) {
    /// special case for cluster.lip6.fr 
    // supply n_max n_min l_max l_min
    int n_max = 100;
    int n_min = 71;
    int l_max = 31;
    int l_min = 31;
    // l = atof(argv[2]);

    printf("n_max=%d,  n_min=%d, l_max=%d, l_min=%d\n",
 	   n_max, n_min, l_max, l_min);
    // variable file name

    char file_name[43];
    snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_%d_%d_%d_stats.txt",
	     n_max, n_min, l_max, l_min);
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
  
  else if (argc == 3){ // ./long_message_attack n l
  // get the input from the user
  n = atoi(argv[1]);
  l = atof(argv[2]);

  // variable file name
  char file_name[36];
  snprintf(file_name, sizeof(file_name), "%s%d_%d_stats.txt", directory_name, (int) n, (int) l);


  /// todo write the value of n and l explicitly in the file 
  //  FILE* fp = fopen("statistics/n_l_stats.txt", "w");
  FILE* fp = fopen(file_name, "w");
  fprintf(fp, "%s", first_line);
	    

  long_message_attack(n, l, fp);

  fclose(fp);

  return 0;
    }

  else if (argc == 2){ // ./long_message_attack n
    // we test all n1 <= n and i <l<= n/2
    // i is determined by the programmer 
    n = atoi(argv[1]);

    // variable file name
    char file_name[36];
    snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_stats.txt", (int) n);


    
    
    FILE* fp = fopen(file_name, "w");
    fprintf(fp, "%s", first_line);
    fclose(fp);

    int l_max = 26; // after this the program consumes more than the available ram
    for (int n1 = 2; n1 < n; ++n1){
      // when l > n/2 then the we expect to be a cylce in phase I
      for (int l=1; (l<=l_max && l<= (n1>>1) ); ++l){
	// printf("n=%d, l=%d\n", n1, l );
	// opening file multiple time to results as soon we as we have it
	FILE* fp = fopen(file_name, "a");
	long_message_attack(n1, l, fp);
	fclose(fp);
	// puts("");
      }
    }
    
    
  } else if (argc == 5) { // ./long_message_attack nmax nmin lmax lmin

    // supply n_max n_min l_max l_min
    int n_max = atoi(argv[1]);
    int n_min = atoi(argv[2]);
    int l_max = atoi(argv[3]);
    int l_min = atoi(argv[4]);
    // l = atof(argv[2]);

    printf("n_max=%d,  n_min=%d, l_max=%d, l_min=%d\n",
 	   n_max, n_min, l_max, l_min);
    // variable file name

    char file_name[43];
    snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_%d_%d_%d_stats.txt",
	     n_max, n_min, l_max, l_min);
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



