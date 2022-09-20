/// we define here how large uint

#ifndef DICT_ELM_TYPE
#define DICT_ELM_TYPE
#include <stdint.h>


/// Add macro that control the output length of the digest

#define UINT uint32_t
// #define UINT uint32_t

typedef union {
  /// Wrap state type in a union for easier handling outside the sha256 function 
  char bytes[32];
  unsigned int state[8]; // as used in sha256 code
  uint64_t values[4];

  } digest;


typedef union {
  /// Wrap state type in a union for easier handling outside the sha256 function 
  char bytes[8];
  unsigned int _uint32[2]; // as used in sha256 cod
  uint64_t _uint64[1];

  } dict_key;

  
#endif




