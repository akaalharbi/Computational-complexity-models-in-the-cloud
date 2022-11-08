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
#include <endian.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include <sys/time.h>
#include <omp.h>

// mask for distinguished point
#define DIST_MASK 0x7 // test if the last three bits are zeros

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///


// was there a cycle in PHASE I
int is_there_duplicate = 0;
int idx_cycle = -1;



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
  size_t nmsg_rec = 0; // number of messages recorded

  while (nmsg_rec < n_of_blocks){
    // number of (distinguished messages == n_of_blocks)
    sha256_process_x86_single(state, M);

    // get the digest from the state
    truncate_state32bit_get_digest(digest, state, n_of_bits);
    
    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///

    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

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
  }
  ///             DONE WITH HASHING THE LONG MESSAG            ///
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
    printf("Thread %02d\n", omp_get_thread_num());


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
      sha256_process_x86_single(state_priv, random_message_priv);

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


  /* if (argc == 1) { */
  /*   /// special case for cluster.lip6.fr  */
  /*   // supply n_max n_min l_max l_min */
  /*   int n_max = 100; */
  /*   int n_min = 71; */
  /*   float l = 33.5; */
  /*   // l = atof(argv[2]); */

  /*   printf("n_max=%d,  n_min=%d, l_max=%f, l_min=%f\n", */
  /* 	   n_max, n_min, l, l); */
  /*   // variable file name */

  /*   char file_name[43]; */
  /*   snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_%d_one_l_stats.txt", */
  /* 	     n_max, n_min); */
  /*   FILE* fp = fopen(file_name, "w"); */
  /*   fprintf(fp, "%s", first_line);    fclose(fp); */
    
  /*   /// loop over n_min <= n <= n_max, l_min <= l <= l_max */
  /*   // opening file multiple time to results as soon we as we have it */
   
  /*   for (int n1=n_min; n1<=n_max; ++n1){ */
  /*     FILE* fp = fopen(file_name, "a"); */
  /*     long_message_attack(n1, l, fp); */
  /*     fclose(fp);       */
  /*   } */


  /* 	// puts(""); */
  /* } */
    

  
  if (argc == 3){ // ./long_message_attack n l
  // get the input from the user
  n = atoi(argv[1]);
  l = atof(argv[2]);

  // variable file name
  char file_name[40];
  snprintf(file_name, sizeof(file_name), "%s%d_%d_stats.txt", directory_name, (int) n, (int) l);
  printf("filename=%s\n", file_name);

  /// todo write the value of n and l explicitly in the file 
  //  FILE* fp = fopen("statistics/n_l_stats.txt", "w");
  FILE* fp = fopen(file_name, "w");
  fprintf(fp, "%s", first_line);
	    

  long_message_attack(n, l, fp);

  fclose(fp);

  return 0;
    }

  /* else if (argc == 2){ // ./long_message_attack n */
  /*   // we test all n1 <= n and i <l<= n/2 */
  /*   // i is determined by the programmer  */
  /*   n = atoi(argv[1]); */

  /*   // variable file name */
  /*   char file_name[36]; */
  /*   snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_stats.txt", (int) n); */


    
    
  /*   FILE* fp = fopen(file_name, "w"); */
  /*   fprintf(fp, "%s", first_line); */
  /*   fclose(fp); */

  /*   int l_max = 26; // after this the program consumes more than the available ram */
  /*   for (int n1 = 2; n1 < n; ++n1){ */
  /*     // when l > n/2 then the we expect to be a cylce in phase I */
  /*     for (int l=1; (l<=l_max && l<= (n1>>1) ); ++l){ */
  /* 	// printf("n=%d, l=%d\n", n1, l ); */
  /* 	// opening file multiple time to results as soon we as we have it */
  /* 	FILE* fp = fopen(file_name, "a"); */
  /* 	long_message_attack(n1, l, fp); */
  /* 	fclose(fp); */
  /* 	// puts(""); */
  /*     } */
  /*   } */
    
    
  /* } else if (argc == 5) { // ./long_message_attack nmax nmin lmax lmin */

  /*   // supply n_max n_min l_max l_min */
  /*   int n_max = atoi(argv[1]); */
  /*   int n_min = atoi(argv[2]); */
  /*   int l_max = atoi(argv[3]); */
  /*   int l_min = atoi(argv[4]); */
  /*   // l = atof(argv[2]); */

  /*   printf("n_max=%d,  n_min=%d, l_max=%d, l_min=%d\n", */
  /* 	   n_max, n_min, l_max, l_min); */
  /*   // variable file name */

  /*   char file_name[43]; */
  /*   snprintf(file_name, sizeof(file_name), "statistics_parallel/%d_%d_%d_%d_stats.txt", */
  /* 	     n_max, n_min, l_max, l_min); */
  /*   FILE* fp = fopen(file_name, "w"); */
  /*   fprintf(fp, "%s", first_line);    fclose(fp); */
    
  /*   /// loop over n_min <= n <= n_max, l_min <= l <= l_max */
  /*   for (int n1=n_min; n1 <= n_max; ++n1){ */
  /*     // when l > n/2 then the we expect to be a cylce in phase I */
  /*     for (int l=l_min; (l<=l_max && l<= (n1>>1) ); ++l){ */
  /* 	// opening file multiple time to results as soon we as we have it */
  /* 	FILE* fp = fopen(file_name, "a"); */
  /* 	long_message_attack(n1, l, fp); */
  /* 	fclose(fp); */
  /* 	// puts(""); */
  /*     } */
  /*   } */


  /* } */
  
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



