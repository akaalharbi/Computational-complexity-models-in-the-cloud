#ifndef SHA256_SINGLE
#define SHA256_SINGLE
#include <stdint.h>
void sha256_single(uint32_t state[8], const uint8_t data[]);
#endif
