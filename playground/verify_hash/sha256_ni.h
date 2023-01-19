#ifndef SHA256_NI_WRAPPER

#define SHA256_NI_WRAPPER
#include <stdint.h>

void sha256_process_x86(uint32_t state[8], const uint8_t data[],
                        uint32_t length);



#endif
