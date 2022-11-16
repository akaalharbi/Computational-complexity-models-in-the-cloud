/// we define here how large uint

#ifndef LONG_MESSAGE_CONFIG
#define LONG_MESSAGE_CONFIG
#include <stdint.h>

// dead weight they have to be reomved
#define UINT uint64_t // this
#define BYTE unsigned char

/// From our experiments this number gives the best tradeoff
/// between memory and cpu time
#define FILLING_RATE 0.9

// depending on avx register length, on my laptop 256
//  Change them together, @todo write a code to automate writing these values
// ---------------------- |
#define AVX_SIZE  256 //  |
#define ALIGNMENT 32  //  |
// ---------------------- |

#endif




