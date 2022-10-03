/// The purpose of this file is to benchmark
/// sha256 calls
/// dictionary adding and probing an element giving a load value in (0, 1)


#include "sha256.h"
#include <stddef.h>
#include <stdio.h>
#include "dict.h"
#include <stdlib.h>
#include <sys/time.h>
#include "shared.h"
#include <math.h>


// was there a cycle in PHASE I
int is_there_duplicate = 0;
int idx_cycle = -1;



float bench_mark_sha256(){
  size_t n_of_blocks = 1<<25;
  
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  // uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  SHA256_CTX ctx;
  sha256_init(&ctx);


  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);
  for (size_t i=0; i<n_of_blocks; ++i){
    sha256_transform(&ctx, M);
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
  printf("elapsed=%fsec, %f hash/sec≈2^%f \n", elapsed, hash_per_sec, log2(hash_per_sec));
  return hash_per_sec;
}


void filling_rate_time(size_t n_of_blocks, float alpha){
  // size_t n_of_blocks = 1<<25;
  // the dictionary has size 2^26
  dict* d = dict_new(n_of_blocks);

  size_t N = (size_t) (n_of_blocks<<1) * alpha;
  printf("N=%lu=2^%f\n", N, log2(N));
  printf("dict has %lu slots = 2^%f slots", d->nslots, log2(d->nslots));

    struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;

  BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);
  // store the hash value in this variable
  uint64_t digest[2] = {0, 0};
  // INIT SHA256 
  SHA256_CTX ctx;
  sha256_init(&ctx);


  // hash a long message (for now it's a series of zeros)
  gettimeofday(&begin, 0);
  for (size_t i=0; i<N; ++i){
    sha256_transform(&ctx, M);
    truncate_state_get_digest(digest, &ctx, 128);
    dict_add_element_to(d, digest, i);  
  }
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;

  
  printf("dictionary filling %lu elements took %fsec \n", N,  elapsed);
  printf("i.e. %f elm/sec≈2^%felm/sec\n", (float) N / elapsed, log2((float) N / elapsed));
  puts("--------------------");
  // dummy variable so the compiler won't optimized the loop below
  size_t values = 0; 
  gettimeofday(&begin, 0);
  for (size_t i=0; i<N; ++i){
    sha256_transform(&ctx, M);
    truncate_state_get_digest(digest, &ctx, 128);
    values += dict_get_value(d, digest);
  }
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  printf("dictionary lookup %lu elements took %fsec \n", N,  elapsed);
  printf("i.e. %f elm/sec≈2^%felm/sec\n", (float) N / elapsed, log2((float) N / elapsed));

  // 
}

 
int main(int argc, char* argv[]){
  bench_mark_sha256();
  filling_rate_time(1<<25, 0.95);
}


// alpha N log2(N) (avg insert_time (elm/sec)) (avg lookup_time (elm/sec))
