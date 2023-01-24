/*******************************************************************************
  Copyright (c) 2022, Intel Corporation

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

      * Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of Intel Corporation nor the names of its contributors
        may be used to endorse or promote products derived from this software
        without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#ifndef IMB_ASM_AVX2_T1_H
#define IMB_ASM_AVX2_T1_H

#include <inttypes.h>
// #include "intel-ipsec-mb.h"
// #include "ipsec_ooo_mgr.h"

// #define AVX2_NUM_SHA256_LANES   8
#define AVX512_NUM_SHA256_LANES 16
#define SHA256_DIGEST_WORD_SIZE   4
#define NUM_SHA_256_DIGEST_WORDS 8

#define DECLARE_ALIGNED(decl, alignval) decl __attribute__((aligned(alignval)))

#define SHA256_DIGEST_SZ (NUM_SHA_256_DIGEST_WORDS * AVX512_NUM_SHA256_LANES)

typedef struct {
        DECLARE_ALIGNED(uint32_t digest[SHA256_DIGEST_SZ], 32);
        const uint8_t *data_ptr[AVX512_NUM_SHA256_LANES];
} SHA256_ARGS;

/* SHA */
void call_sha256_oct_avx2_from_c(SHA256_ARGS *args, uint32_t size_in_blocks);

#endif /* IMB_ASM_AVX2_T1_H */
