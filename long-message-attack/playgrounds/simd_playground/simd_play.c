// __m256i _mm256_load_si256 (__m256i const * mem_addr)
#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#define KEY_TYPE uint64_t
#define ALIGNMENT 32
#define NVECT 2
#define STEP (NVECT*ALIGNMENT/sizeof(KEY_TYPE)) // for the while loop

void print_array(uint64_t* a, int nelements){
  for (int i = 0; i<nelements; ++i)
    printf("%lx, ", a[i]);
  puts("");
}

int main(int argc, char* argv[]){
  int nelements = 8;
  printf("step=%lu\n", STEP);
  uint64_t* A = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t)));
  uint64_t* B = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t)));
  uint64_t* C = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t)));
  for (int i = 0; i<nelements; ++i){
    A[i] = i;
    B[i] = 2*i;
  }

  __m256i a =  _mm256_stream_load_si256((__m256i*) A);
  __m256i b =  _mm256_stream_load_si256((__m256i*) B);
  
  __m256i c = _mm256_add_epi64(a, b);
  _mm256_store_si256((void*)C, c);

  print_array(C, nelements);

  
  
  
}
