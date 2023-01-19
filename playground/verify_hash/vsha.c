#include "vsha256.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#define SIMD_SIZE 8

void vsha256_init(uint32_t **states)
{
  *states = (uint32_t*) malloc(sizeof(uint32_t)*8*SIMD_SIZE);
  
};

void vsha256(uint32_t state[8], const uint8_t data[],
                        uint32_t length);
