
#ifndef SHA256_INTEL_AVX512
#define SHA256_INTEL_AVX512
#include <stdint.h>
#include "arch_avx512_type1.h"

uint32_t *sha256_multiple_x16(uint8_t msg[16][64], uint32_t state[8]);
void call_sha256_x16_avx512_from_c(SHA256_ARGS *args, uint32_t size_in_blocks);

#endif
