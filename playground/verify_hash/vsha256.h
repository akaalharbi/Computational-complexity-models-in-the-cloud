#ifndef SHA256_NI_WRAPPER
#define SHA256_NI_WRAPPER
#include <stdint.h>



void vsha256_init(uint32_t **states);

void vsha256(uint32_t state[8], const uint8_t data[],
                        uint32_t length);


#endif
