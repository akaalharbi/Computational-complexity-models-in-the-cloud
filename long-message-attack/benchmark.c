/// The purpose of this file is to benchmark
/// sha256 calls
/// dictionary adding and probing an element giving a load value in (0, 1)
/// todo generic code that works with predefined FILLING_RATE
/// the current code works with FILLING_RATE = 0.5


#include "sha256-x86.h"
#include <stddef.h>
#include <stdio.h>
#include "dict.h"
#include <stdlib.h>
#include <sys/time.h>
#include "shared.h"
#include <math.h>
#include <omp.h>

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



void filling_rate_time(size_t n_of_blocks, float alpha, FILE* fp){
  // size_t n_of_blocks = 1<<25;
  // the dictionary has size 2^26
  dict* d = dict_new(n_of_blocks);

  size_t N = (size_t) (n_of_blocks<<1) * alpha;
  /* printf("N=%lu=2^%f\n", N, log2(N)); */
  /* printf("dict has %lu slots = 2^%f slots", d->nslots, log2(d->nslots)); */
  fprintf(fp, "%.2f, ", alpha);
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;
  double elapsed_total=0;
  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  uint64_t digest[2] = {0, 0};
  // INIT SHA256
  uint32_t state[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };


  // --------------------------------------///
  //// INSERTION TIMING FINE TUNED
  for (size_t i=0; i<N; ++i){
    sha256_process_x86_single(state, M);
    truncate_state32bit_get_digest(digest, state, 128);

    gettimeofday(&begin, 0);
    dict_add_element_to(d, digest, i);
    gettimeofday(&end, 0);    
    seconds = end.tv_sec - begin.tv_sec;
    microseconds = end.tv_usec - begin.tv_usec;
    elapsed = seconds + microseconds*1e-6;
    elapsed_total += elapsed;

  }

  // how wasteful !
  elapsed = elapsed_total;
  elapsed_total = 0;
  
  printf("dictionary filling %lu elements took %fsec \n", N,  elapsed);
  printf("i.e. %f elm/sec≈2^%felm/sec\n", (float) N / elapsed, log2((float) N / elapsed));
  fprintf(fp, "%felm/sec, %fprobes/elm, ", (float) N / elapsed, ((float)d->nprobes_insert)/N);
  // edit base_alpha only in the first call of the function
  if (base_alpha == 0)
    base_alpha = alpha;
 


  /// LOOKUP TIMING FINE TUNED
  size_t values = 0; 

  for (size_t i=0; i<N; ++i){
    // random numbers
    sha256_process_x86_single(state, M);

    truncate_state32bit_get_digest(digest, state, 128);

    gettimeofday(&begin, 0);
    values += dict_get_value(d, digest);
    gettimeofday(&end, 0);
    seconds = end.tv_sec - begin.tv_sec;
    microseconds = end.tv_usec - begin.tv_usec;
    elapsed = seconds + microseconds*1e-6;
    elapsed_total += elapsed;
  }

  // how wasteful !
  elapsed = elapsed_total;
  elapsed_total = 0;

  printf("dictionary lookup %lu elements took %fsec \n", N,  elapsed);
  // fill base_lookup_rate only in the first call
  if (base_lookup_rate==0)
    base_lookup_rate = ((float) N) / elapsed;

  float new_lookup_rate = ((float) N) / elapsed;
  // how many bits have we gain in memory
  float gain = log2(alpha/base_alpha);
  // how many bits have we lost in performance
  float loss = log2(new_lookup_rate/base_lookup_rate);
  printf("i.e. %f elm/sec≈2^%felm/sec\n", new_lookup_rate, log2((float) new_lookup_rate));
  fprintf(fp, "%felm/sec, %fprobes/elm, %fbits\n", new_lookup_rate, ((float) d->nprobes_lookup)/N, gain+loss);
  puts("--------end---------");
  dict_free(d);
  free(d);
}

 
int main(int argc, char* argv[]){
  /// Planning
  /// open file named dict_benchmark in log
  

  benchmark_sha256_x86();

  // benchmark parallel sha256
  benchmark_sha256_x86_parallel();

  
  /* // benchmark filling rate */
  /* size_t nelements = 1<<25; */
  /* FILE* fp = fopen("log/benchmark_dict", "w"); */
  /* fprintf(fp, "alpha, insert, nprobes_insert,  lookup, nprobes_lookup, total_gain\n" */
  /* 	  "N=%lu\n", nelements); */
  /* fclose(fp); */
  
  /* for (float i=0.5; i<0.99; i += 0.01){ */
  /*   FILE* fp = fopen("log/benchmark_dict", "a"); */
  /*   filling_rate_time(nelements, i, fp); */
  /*   fclose(fp); */
  /* } */
}


// alpha N log2(N) (avg insert_time (elm/sec)) (avg lookup_time (elm/sec))
