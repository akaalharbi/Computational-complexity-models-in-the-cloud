/// The purpose of this file is to benchmark
/// sha256 calls
/// dictionary adding and probing an element giving a load value in (0, 1)
/// todo generic code that works with predefined FILLING_RATE
/// the current code works with FILLING_RATE = 0.5


#include "sha256-x86.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
//#include "dict.h"
#include <stdlib.h>
#include <sys/time.h>
#include "shared.h"
#include <math.h>
#include <omp.h>
#include "vsha256.h"
#define BYTE unsigned char

// was there a cycle in PHASE I
int is_there_duplicate = 0;
int idx_cycle = -1;

float base_alpha = 0;
float base_lookup_rate = 0;



float benchmark_sha256_x86(){

  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  
  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };



  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_process_x86_single(state, M);
    // truncate_state_get_digest(digest, &ctx, n_of_bits);

    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///
    /// time. 
    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

  }



  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) n_of_blocks / elapsed;
  printf("sha256-x86\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}


float benchmark_sha256_x86_parallel(){

  
  size_t n_of_blocks = 1<<25;
  int nthreads = omp_get_max_threads();
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  #pragma omp parallel
  {
  BYTE M[64]  __attribute__ ((aligned (32)))  = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  
  uint32_t state[8] __attribute__ ((aligned (32))) = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  } ;



  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);  
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_process_x86_single(state, M);
     }
  }

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) (nthreads*n_of_blocks) / elapsed;
  printf("parallel sha256-x86\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}

float benchmark_sha256_x86_modified(){

  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  
  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };



  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_process_x86_single_modified(state, M);
    // truncate_state_get_digest(digest, &ctx, n_of_bits);

    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///
    /// time. 
    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

  }



  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) n_of_blocks / elapsed;
  printf("sha256-x86-modified\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}

float benchmark_sha256_x86_modified_parallel(){

  
  size_t n_of_blocks = 1<<25;
  int nthreads = omp_get_max_threads();
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  #pragma omp parallel
  {
  BYTE M[64] __attribute__ ((aligned (32)))  = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  
  uint32_t state[8]  __attribute__ ((aligned (32)))   = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };



  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);  
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_process_x86_single_modified(state, M);
     }
  }

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) (nthreads*n_of_blocks) / elapsed;
  printf("parallel sha256-x86-modified\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}

// sha256 vecotr implementation
float benchmark_vsha256(){

  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  // store 8x256-bits hashes
  u32 AVX_ALIGNED h[8][8];
  vsha256_init(h); // 
  u32 AVX_ALIGNED msg[16][8];
  for (int i = 0; i < 8; i++)
    msg[0][i] = 0x00010203 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[1][i] = 0x04050607 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[2][i] = 0x08090a0b + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[3][i] = 0x0c0d0e0f + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[4][i] = 0x10111213 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[5][i] = 0x14151617 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[6][i] = 0x18191a1b + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[7][i] = 0x1c1d1e1f + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[8][i] = 0x20212223 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[9][i] = 0x24252627 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[10][i] = 0x28292a2b + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[11][i] = 0x2c2d2e2f + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[12][i] = 0x30313233 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[13][i] = 0x34353637 + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[14][i] = 0x38393a3b + 0x10101010 * i;
  for (int i = 0; i < 8; i++)
    msg[15][i] = 0x3c3d3e3f + 0x10101010 * i;


  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  




  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);
  for (size_t i=0; i<(n_of_blocks>>3); ++i){
    vsha256_transform(h, msg);
    // truncate_state_get_digest(digest, &ctx, n_of_bits);

    
    /// ------------ DISTINGUISHED POINTS ------------------------- ///
    /// If distinguished points feature was enabled  during compile ///
    /// time. 
    #ifdef DISTINGUISHED_POINTS
    // we skip hashes
    if ( (digest[0]&DIST_MASK) != 0) 
      continue; // skip this element
    #endif

  }



  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) n_of_blocks / elapsed;
  printf("vsha256\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}


// sha256 vecotr implementation
float benchmark_vsha256_parallel(){

  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;


  #pragma omp parallel
  {

    // store 8x256-bits hashes
    u32 AVX_ALIGNED h[8][8];
    vsha256_init(h); // 
    u32 AVX_ALIGNED msg[16][8];
    for (int i = 0; i < 8; i++)
      msg[0][i] = 0x00010203 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[1][i] = 0x04050607 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[2][i] = 0x08090a0b + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[3][i] = 0x0c0d0e0f + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[4][i] = 0x10111213 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[5][i] = 0x14151617 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[6][i] = 0x18191a1b + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[7][i] = 0x1c1d1e1f + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[8][i] = 0x20212223 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[9][i] = 0x24252627 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[10][i] = 0x28292a2b + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[11][i] = 0x2c2d2e2f + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[12][i] = 0x30313233 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[13][i] = 0x34353637 + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[14][i] = 0x38393a3b + 0x10101010 * i;
    for (int i = 0; i < 8; i++)
      msg[15][i] = 0x3c3d3e3f + 0x10101010 * i;


    // store the hash value in this variable
    // uint64_t digest[2] = {0, 0};
    // INIT SHA256 
  




    // hash a long message (for now it's a series of zeros)
    gettimeofday(&begin, 0);
    for (size_t i=0; i<(n_of_blocks>>3); ++i){
      vsha256_transform(h, msg);
      // truncate_state_get_digest(digest, &ctx, n_of_bits);

    
      /// ------------ DISTINGUISHED POINTS ------------------------- ///
      /// If distinguished points feature was enabled  during compile ///
      /// time. 
#ifdef DISTINGUISHED_POINTS
      // we skip hashes
      if ( (digest[0]&DIST_MASK) != 0) 
	continue; // skip this element
#endif

    }

  }

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  int nthreads = omp_get_max_threads();
  float hash_per_sec = (float) nthreads*n_of_blocks / elapsed;
  printf("vsha256_parallel\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}






/* void filling_rate_time(size_t n_of_blocks, float alpha, FILE* fp){ */
/*   // size_t n_of_blocks = 1<<25; */
/*   // the dictionary has size 2^26 */
/*   dict* d = dict_new(n_of_blocks); */

/*   size_t N = (size_t) (n_of_blocks<<1) * alpha; */
/*   printf("N=%lu=2^%f\n", N, log2(N)); */
/*   printf("dict has %lu slots = 2^%f slots", d->nslots, log2(d->nslots)); */
/*   fprintf(fp, "%.2f, ", alpha); */
/*   struct timeval begin, end; */
/*   long seconds = 0; */
/*   long microseconds = 0; */
/*   double elapsed = 0; */
/*   double elapsed_total=0; */
/*   BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512); */
/*   // store the hash value in this variable */
/*   uint64_t digest[2] = {0, 0}; */
/*   // INIT SHA256 */
/*   uint32_t state[8] = { */
/*     0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, */
/*     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 */
/*   }; */


/*   // --------------------------------------/// */
/*   //// INSERTION TIMING FINE TUNED */
/*   for (size_t i=0; i<N; ++i){ */
/*     sha256_process_x86_single(state, M); */
/*     truncate_state32bit_get_digest(digest, state, 128); */

/*     gettimeofday(&begin, 0); */
/*     dict_add_element_to(d, digest); */
/*     gettimeofday(&end, 0);     */
/*     seconds = end.tv_sec - begin.tv_sec; */
/*     microseconds = end.tv_usec - begin.tv_usec; */
/*     elapsed = seconds + microseconds*1e-6; */
/*     elapsed_total += elapsed; */

/*   } */

/*   // how wasteful ! */
/*   elapsed = elapsed_total; */
/*   elapsed_total = 0; */
  
/*   printf("dictionary filling %lu elements took %fsec \n", N,  elapsed); */
/*   printf("i.e. %f elm/sec≈2^%felm/sec\n", (float) N / elapsed, log2((float) N / elapsed)); */
/*   fprintf(fp, "%felm/sec, %fprobes/elm, ", (float) N / elapsed, ((float)d->nprobes_insert)/N); */
/*   // edit base_alpha only in the first call of the function */
/*   if (base_alpha == 0) */
/*     base_alpha = alpha; */
 


/*   /// LOOKUP TIMING FINE TUNED */
/*   // size_t values = 0;  */
/*   #define NSIMD_SHA 4 */
/*   // use simd to create 8 hashes simultanously */
/*   uint32_t state_init_priv[8] = { */
/*     0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, */
/*     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 */
/*   }; */
 
/*   BYTE random_message_priv[NSIMD_SHA][64] = {0}; */
/*   uint64_t digest_priv[NSIMD_SHA][2] = {0}; */
/*   uint64_t lookup_keys_priv[NSIMD_SHA] = {0}; */
/*   uint32_t state_priv[NSIMD_SHA][8] = {0}; */
/*   size_t found_keys_priv[NSIMD_SHA] = {0}; */

  
/*   for (size_t i=0; i<(N/NSIMD_SHA); ++i){ */
/*     // random numbers */
/*     for (int i=0; i<NSIMD_SHA; ++i) { */
/*       fill_radom_byte_array_get_random(random_message_priv[i], 64); */
/*       // clean previously used values */
/*       memcpy(state_priv[i], state_init_priv, sizeof(state_init_priv)); */
/*       sha256_process_x86_single(state_priv[i], random_message_priv[i]);	 */
/*       truncate_state32bit_get_digest(digest_priv[i], state_priv[i], 128); */
/*       lookup_keys_priv[i] = digest_priv[i][0]; */
/*     } */
/*     gettimeofday(&begin, 0); */
/*     dict_get_values_simd(d, lookup_keys_priv, found_keys_priv); */
/*     gettimeofday(&end, 0); */

/*     seconds = end.tv_sec - begin.tv_sec; */
/*     microseconds = end.tv_usec - begin.tv_usec; */
/*     elapsed = seconds + microseconds*1e-6; */
/*     elapsed_total += elapsed; */
/*   } */

/*   // how wasteful ! */
/*   elapsed = elapsed_total; */
/*   elapsed_total = 0; */

/*   printf("dictionary lookup %lu elements took %fsec \n", N,  elapsed); */
/*   printf("dictionary successful lookups %lu \n", d->nelments_succ_lookup); */
/*   // fill base_lookup_rate only in the first call */
/*   if (base_lookup_rate==0) */
/*     base_lookup_rate = ((float) N) / elapsed; */

/*   float new_lookup_rate = ((float) N) / elapsed; */
/*   // how many bits have we gain in memory */
/*   float gain = log2(alpha/base_alpha); */
/*   // how many bits have we lost in performance */
/*   float loss = log2(new_lookup_rate/base_lookup_rate); */
/*   printf("i.e. %f elm/sec≈2^%felm/sec\n", new_lookup_rate, log2((float) new_lookup_rate)); */
/*   fprintf(fp, "%felm/sec, %fprobes/elm, %fbits\n", new_lookup_rate, ((float) d->nprobes_lookup)/N, gain+loss); */
/*   puts("--------end---------"); */
/*   dict_free(d); */
/*   free(d); */
/* } */

 
int main(int argc, char* argv[]){
  /// Planning
  /// open file named dict_benchmark in log
  

  benchmark_sha256_x86();

  // benchmark parallel sha256
  benchmark_sha256_x86_parallel();

  benchmark_sha256_x86_modified();
  benchmark_sha256_x86_modified_parallel();
  benchmark_vsha256();
  benchmark_vsha256_parallel();

  /* // benchmark filling rate */
  /* size_t nelements = 1<<25; */

  /* FILE* fp = fopen("log/benchmark_dict", "w"); */
  /* fprintf(fp, "alpha, insert, nprobes_insert,  lookup, nprobes_lookup, total_gain\n" */
  /* 	  "N=%lu\n", nelements); */
  /* fclose(fp); */
  /* //filling_rate_time(nelements, 0.9, fp);   */
  /* for (float i=0.5; i<0.99; i += 0.01){ */
  /*   FILE* fp = fopen("log/benchmark_dict", "a"); */
  /*   filling_rate_time(nelements, i, fp); */
  /*   fclose(fp); */
  /* } */
}


// alpha N log2(N) (avg insert_time (elm/sec)) (avg lookup_time (elm/sec))
