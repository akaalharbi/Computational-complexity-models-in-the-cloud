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

/* void print_m25i(__m256i a, char* text){ */
/*   uint32_t A[8] = {0}; */
/*   _mm256_storeu_si256((__m256i*)A, a); */
/*   printf("%s = ", text); */
/*   for (int i = 0; i<8; ++i) { */
/*     printf("%02x, ", A[i]); */
/*   } */
/*   puts(""); */
/* } */

void print_m25i(__m256i a, char* text){
  uint64_t A[4] = {0};
  _mm256_storeu_si256((__m256i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<4; ++i) {
    printf("%016lx, ", A[i]);
  }
  puts("");
}

void print_m128i(__m128i a, char* text){
  uint32_t A[4] = {0};
  _mm_storeu_si128((__m128i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<4; ++i) {
    printf("%02x, ", A[i]);
  }
  puts("");
}


int main(int argc, char* argv[]){
  /* int nelements = 8; */
  /* printf("step=%lu\n", STEP); */
  /* uint64_t* A = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t))); */
  /* uint64_t* B = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t))); */
  /* uint64_t* C = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t))); */
  /* for (int i = 0; i<nelements; ++i){ */
  /*   A[i] = i; */
  /*   B[i] = 2*i; */
  /* } */

  /* __m256i a =  _mm256_stream_load_si256((__m256i*) A); */
  /* __m256i b =  _mm256_stream_load_si256((__m256i*) B); */
  
  /* __m256i c = _mm256_add_epi64(a, b); */
  /* _mm256_store_si256((void*)C, c); */

  /* print_array(C, nelements); */

  /* // test gather */
  /* puts("GATHER TEST: I give up on using gather x)"); */
  /* nelements = 9; // we write only on the first 4, the others are for testing */
  /* uint64_t data[9] ={0x29, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19}; */
  /* uint64_t* loaded = (uint64_t*) aligned_alloc(ALIGNMENT, nelements*(sizeof(uint64_t))); */

  /* __m256i indices = _mm256_set_epi64x(3, */
  /* 				      2, */
  /* 				      1, */
  /* 				      0); */
  /* __m256i array = _mm256_i64gather_epi64(data, indices, 2);  */
  /* puts("hi");  */

  //_mm256_store_si256((void*)loaded, array);
  /* _mm256_store_si256((void*)loaded, array); */
  /* puts("loaded:"); */
  /* print_array(loaded, nelements); */


  ///////////////////////////////////////////////////////////////////////////
  /// TEST alignr from 32 to 64 bit




  /* STATE0 = _mm256_alignr_epi8(TMP, STATE1, 8);    /\* ABEF *\/ */
  /* STATE1 = _mm_blend_epi16(STATE1, TMP, 0xF0); /\* CDGH *\/ */


  
  /* uint32_t A[16] __attribute__((aligned(32))) = {0x1a, 0xb, 0xc, 0xd, */
  /* 						 0xa, 0xb, 0xc, 0xd, */
  /* 						 0xe, 0xf, 0x10, 0x11, */
  /* 						 0xe, 0xf, 0x10, 0x11}; */


  /*  __m256i STATE0, STATE1; */
  /* __m256i TMP = _mm256_load_si256( (__m256i*) &A[0]); */

  /* STATE1 = _mm256_load_si256( (__m256i*) &A[8]); */
  /* print_m25i(STATE1, "STATE1");   */
  /* print_m25i(TMP,    "TMP   "); */
  
  /* TMP = _mm256_shuffle_epi32(TMP, 0xB1);          /\* CDAB *\/   */
  /* STATE1 = _mm256_shuffle_epi32(STATE1, 0x1B);    /\* EFGH *\/   */

  /* print_m25i(TMP,    "shuffle TMP   "); */
  /* print_m25i(STATE1, "shuffle STATE1"); */
  
  /* STATE0 = _mm256_alignr_epi8(TMP, STATE1, 8);    /\* ABEF *\/ */
  /* STATE1 = _mm256_blend_epi16(STATE1, TMP, 0xF0); /\* CDGH *\/ */


  /* print_m25i(STATE0, "alignr STATE0"); */
  /* print_m25i(STATE1, "blend  STATE1"); */
  /* puts("G-> 10, H->11"); */


  /* uint32_t B[8] __attribute__((aligned(32))) = {0xa, 0xb, 0xc, 0xd, 0xe, 0xf, 0x10, 0x11}; */
  /* __m128i tmp = _mm256_castsi256_si128( TMP); //_mm_load_si128( (__m128i*) &B[0]); */
  /* print_m128i(tmp,    "tmp"); */
  /* tmp = _mm_shuffle_epi32(tmp, 0xB1);          /\* CDAB *\/   */
  /* print_m128i(tmp,    "tmp after"); */


  
  uint64_t AMIZERO[4] = {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
  __m256i amizero = _mm256_loadu_si256((__m256i*) AMIZERO);
  print_m25i(amizero, "amizero begin");
  int is_it = _mm256_testz_si256(amizero, amizero);
  printf("is_it=%d\n", is_it);
  #ifdef __SIZEOF_INT128__
  puts("we have 128bit support");
  #endif 
}















































