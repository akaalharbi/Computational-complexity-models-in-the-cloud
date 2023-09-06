/// All sha256 implmentations used in long-message-attack
/// must confirm to this uniform interface

#ifndef SHA256_UNIFIED_INTERFACE
#define SHA256_UNIFIED_INTERFACE



#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#if defined(__GNUC__)
# include <stdint.h>
# include <x86intrin.h>
#endif

/* Microsoft supports Intel SHA ACLE extensions as of Visual Studio 2015 */
#if defined(_MSC_VER)
# include <immintrin.h>
# define WIN32_LEAN_AND_MEAN
# include <Windows.h>
typedef UINT32 uint32_t;
typedef UINT8 uint8_t;
#endif

void hash_single(WORD_TYPE state[8], const u8 data[]);

#define HASH_INIT_STATE  0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,\
                         0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19

#endif
