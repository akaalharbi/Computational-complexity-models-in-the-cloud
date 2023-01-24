
#ifndef SHA256_INTEL_AVX2
#define SHA256_INTEL_AVX2

#include <stdint.h>


uint32_t *sha256_multiple_oct(uint8_t msg[16][64]);

#endif
