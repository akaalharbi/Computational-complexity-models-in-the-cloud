
#ifndef SHA256_INTEL_AVX2
#define SHA256_INTEL_AVX2
#include <stdint.h>



uint32_t *sha256_multiple_oct(uint8_t msg[16][64], uint32_t const state[8]);


#endif


#ifndef SHA256_INTEL_AVX512
#define SHA256_INTEL_AVX512


#include "arch_avx512_type1.h"
uint32_t *sha256_multiple_x16(uint8_t msg[16][64]);
void sha256_mb_init_digest_avx512(uint32_t *digest);


#endif
