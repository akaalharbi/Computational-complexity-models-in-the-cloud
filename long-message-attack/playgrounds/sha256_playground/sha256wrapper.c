#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

// for sha256_mult_avx, spoiler alert it doesn't work
#include "intel-headers/ipsec_ooo_mgr.h"
#include "sha256-x86.h"
#include "intel-headers/arch_avx_type1.h"


// for the benchmark
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>


#define SHA256_LEN 32
#define SHA256_INP_LEN 64

extern void sha256_block_avx(const void *, void *);


// ------------- UTIL --------------------- //



void print_char_array(unsigned char* A, size_t len){
  for (size_t i=0; i<len; i++) {
    printf("%02x", A[i]);
  }
  puts("");
}

float benchmark_sha256_x86(){

  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  unsigned char M[64] = {0}; // long_message_zeros(n_of_blocks*512);
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
    sha256_block_avx(state, M);
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
  printf("sha256-intel\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
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
  unsigned char M[64] = {0}; // long_message_zeros(n_of_blocks*512);
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
    sha256_block_avx(state, M);
     }
  }

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  float hash_per_sec = (float) (nthreads*n_of_blocks) / elapsed;
  printf("parallel sha256-intel\nelapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}



int main(int argc, char* argv[]){
  unsigned char data[64] = {0};
  unsigned char state[SHA256_LEN] = {0};

  puts("zeros digest");
  
  sha256_block_avx(data, state);
  print_char_array(state, SHA256_LEN);





  //data[64] = {0};
  //state[32] = {0};
  memset(data, 0, sizeof(data));
  memset(state, 0, sizeof(state));
  // puts("check it is set to zero");
  // print_char_array(state, SHA256_LEN);
  
  //sha256_block_avx(data, state);
  sha256_process_x86_single((uint32_t*)state, data);
  print_char_array(state, SHA256_LEN);

  puts("Benchmark");
  benchmark_sha256_x86();
  benchmark_sha256_x86_parallel();
  
  puts("--------");
  unsigned char mul_state[16*SHA256_LEN] = {0};
  unsigned char mul_data[16*SHA256_INP_LEN] = {0};


  SHA256_ARGS mul_state_data = {0};
  // SHA256_ARGS:
  //   UINT128 digest[8];  // transposed digests
  //   UINT8  *data_ptr[4];
  //

  // void sha_256_mult_avx(SHA256_ARGS *args, UINT64 num_blocks);
  // arg 1 : STATE    : pointer args
  // arg 2 : INP_SIZE : size of data in blocks (assumed >= 1)
  //call_sha_256_mult_avx_from_c(&mul_state_data, 1);

  // print_char_array((unsigned char*) mul_state_data.digest, 16*SHA256_LEN);

}
