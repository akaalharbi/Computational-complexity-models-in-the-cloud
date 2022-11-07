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
    printf("0x%lx, ", a[i]);
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

  // test gather
  puts("GATHER TEST: I give up on using gather x)");
  nelements = 9; // we write only on the first 4, the others are for testing
  uint64_t data[9] ={0x29, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19};
  uint64_t* loaded = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t)));

  __m256i indices = _mm256_set_epi64x(3,
				      2,
				      1,
				      0);
  __m256i array = _mm256_i64gather_epi64(data, indices, 2); 
  puts("hi"); 

  //_mm256_store_si256((void*)loaded, array);
  _mm256_store_si256((void*)loaded, array);
  puts("loaded:");
  print_array(loaded, nelements);

  
}
