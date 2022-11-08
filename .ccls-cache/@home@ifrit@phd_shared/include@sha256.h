/// All sha256 implmentations used in long-message-attack
/// must confirm to this uniform interface

#ifndef SHA256_UNIFIED_INTERFACE
#define SHA256_UNIFIED_INTERFACE



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

void sha256_process_x86_single(uint32_t state[8], const uint8_t data[]);

#endif
