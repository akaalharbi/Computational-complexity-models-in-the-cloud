/// we define here how large uint

#ifndef UINT_TYPE
#define UINT_TYPE
#include <stdint.h>


/// Add macro that control the output length of the digest

#define UINT uint64_t
// #define UINT uint32_t

typedef union {
  /// Wrap state type in a union for easier handling outside the sha256 function 
  char bytes[64];
  uint32_t state[8]; // as used in sha256 code
  uint64_t values[4];

  } digest;

#endif
