/// we define here how large uint

#ifndef LONG_MESSAGE_CONFIG
#define LONG_MESSAGE_CONFIG
#include <stdint.h>

// #define UINT uint32_t 
#define UINT uint64_t // this
#define BYTE unsigned char
/// From our experiments this number gives the best tradeoff
/// between memory and cpu time
#define FILLING_RATE 0.9
#define SIMD_LEN 4


#endif




