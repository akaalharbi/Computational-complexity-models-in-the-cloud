// fast sha256 implementation on intel processeors that supports sha_ni
// source: https://github.com/noloader/SHA-Intrinsics/blob/master/sha256-x86.c
// 

#include "sha256.h"
#include <emmintrin.h>

#define m256i_TO_2m128i(VEC, HI, LO) HI = _mm256_castsi256_si128(VEC);\
  LO = _mm256_extracti128_si256(VEC, 1);\
  


#define two_m128i_to_m256(HI, LO, VEC) VEC = _mm256_castsi128_si256 (LO);	\
  VEC = _mm256_insertf128_si256 (VEC, HI, 1);\
  

 

/* Process multiple blocks. The caller is responsible for setting the initial */
/*  state, and the caller is responsible for padding the final block.        */

void sha256_2way(uint32_t state[16], const uint8_t data[128])
{
  //////////////////////////////////////////////////////////////////////+
  // Different organization:                                            |
  // Say we wish to Hash M := m0 || m1 || m2 || m3 where |mi| = 128-bit |
  // and L = l0 || l1 || l2 || l3                                       |
  // Then data is expected to be inputted interleaved i.e.              |
  // in the form: m0 || ... || m3 || l0 || ... || l3                    |
  //           || m4 || ... || m7 || l4 || ... || l7                    |
  // Also, the state is expected to be interleaved as well              |
  // (128-bit) m-state || (128-bit) l-state || and so on                |
  // In this format, replacing __mm128i (or its variants) to __mm256i ( |
  // or its variant) will work without any modification for most of the |
  // instructions.                                                      |
  //////////////////////////////////////////////////////////////////////+


  /* Casting _mm256i to 2 _mm128i
    lo:  __m128i _mm256_castsi256_si128 (__m256i a)
    hi:  __m128i _mm256_extracti128_si256(data, 1);

     casting 2 __mm128i to__m256i
     
    lo: __m256i _mm256_castsi128_si256 (__m128i a)
    hi: __m256i _mm256_insertf128_si256 (__m256i a, __m128i b, int imm8)

   */
  __m256i STATE0, STATE1;
  __m256i MSG, TMP;
  __m256i MSG0, MSG1, MSG2, MSG3;
  __m256i ABEF_SAVE, CDGH_SAVE;
  const __m256i MASK = _mm256_set_epi64x(0x0c0d0e0f08090a0bULL, 0x0405060700010203ULL,
					 0x0c0d0e0f08090a0bULL, 0x0405060700010203ULL);
  // store first 128-bit of M and L
  __m128i M128_0, M128_1, L128_0, L128_1; 
  /* Load initial values */
  TMP = _mm256_load_si256((const __m256i*) &state[0]);
  STATE1 = _mm256_load_si256((const __m256i*) &state[8]);


  //+ when calling sha_ni, cast 256-bit to 128-bit 
  //
    

  TMP = _mm256_shuffle_epi32(TMP, 0xB1);          /* CDAB */
  STATE1 = _mm256_shuffle_epi32(STATE1, 0x1B);    /* EFGH */
  STATE0 = _mm256_alignr_epi8(TMP, STATE1, 8);    /* ABEF */
  STATE1 = _mm256_blend_epi16(STATE1, TMP, 0xF0); /* CDGH */

  
                                                                                                 
  {
    /* Save current state */
    ABEF_SAVE = STATE0;
    CDGH_SAVE = STATE1;

    /* Rounds 0-3 */
    MSG = _mm256_load_si256((const __m256i*) (data+0));
    MSG0 = _mm256_shuffle_epi8(MSG, MASK);
    MSG = _mm256_add_epi32(MSG0,
			   _mm256_set_epi64x(0xE9B5DBA5B5C0FBCFULL, 0x71374491428A2F98ULL,
					     0xE9B5DBA5B5C0FBCFULL, 0x71374491428A2F98ULL)
			   );
    
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

    /* Rounds 4-7 */
    MSG1 = _mm256_load_si256((const __m256i*) (data+16));
    MSG1 = _mm256_shuffle_epi8(MSG1, MASK);
    MSG = _mm256_add_epi32(MSG1, _mm256_set_epi64x(0xAB1C5ED5923F82A4ULL, 0x59F111F13956C25BULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

    /* Rounds 8-11 */
    MSG2 = _mm256_load_si256((const __m256i*) (data+32));
    MSG2 = _mm256_shuffle_epi8(MSG2, MASK);
    MSG = _mm256_add_epi32(MSG2, _mm256_set_epi64x(0x550C7DC3243185BEULL, 0x12835B01D807AA98ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

    /* Rounds 12-15 */
    MSG3 = _mm256_load_si256((const __m256i*) (data+48));
    MSG3 = _mm256_shuffle_epi8(MSG3, MASK);
    MSG = _mm256_add_epi32(MSG3, _mm256_set_epi64x(0xC19BF1749BDC06A7ULL, 0x80DEB1FE72BE5D74ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG3, MSG2, 4);
    MSG0 = _mm256_add_epi32(MSG0, TMP);
    MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

    /* Rounds 16-19 */
    MSG = _mm256_add_epi32(MSG0, _mm256_set_epi64x(0x240CA1CC0FC19DC6ULL, 0xEFBE4786E49B69C1ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG0, MSG3, 4);
    MSG1 = _mm256_add_epi32(MSG1, TMP);
    MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

    /* Rounds 20-23 */
    MSG = _mm256_add_epi32(MSG1, _mm256_set_epi64x(0x76F988DA5CB0A9DCULL, 0x4A7484AA2DE92C6FULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG1, MSG0, 4);
    MSG2 = _mm256_add_epi32(MSG2, TMP);
    MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

    /* Rounds 24-27 */
    MSG = _mm256_add_epi32(MSG2, _mm256_set_epi64x(0xBF597FC7B00327C8ULL, 0xA831C66D983E5152ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG2, MSG1, 4);
    MSG3 = _mm256_add_epi32(MSG3, TMP);
    MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

    /* Rounds 28-31 */
    MSG = _mm256_add_epi32(MSG3, _mm256_set_epi64x(0x1429296706CA6351ULL,  0xD5A79147C6E00BF3ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG3, MSG2, 4);
    MSG0 = _mm256_add_epi32(MSG0, TMP);
    MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

    /* Rounds 32-35 */
    MSG = _mm256_add_epi32(MSG0, _mm256_set_epi64x(0x53380D134D2C6DFCULL, 0x2E1B213827B70A85ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG0, MSG3, 4);
    MSG1 = _mm256_add_epi32(MSG1, TMP);
    MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

    /* Rounds 36-39 */
    MSG = _mm256_add_epi32(MSG1, _mm256_set_epi64x(0x92722C8581C2C92EULL, 0x766A0ABB650A7354ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG1, MSG0, 4);
    MSG2 = _mm256_add_epi32(MSG2, TMP);
    MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

    /* Rounds 40-43 */
    MSG = _mm256_add_epi32(MSG2, _mm256_set_epi64x(0xC76C51A3C24B8B70ULL, 0xA81A664BA2BFE8A1ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG2, MSG1, 4);
    MSG3 = _mm256_add_epi32(MSG3, TMP);
    MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

    /* Rounds 44-47 */
    MSG = _mm256_add_epi32(MSG3, _mm256_set_epi64x(0x106AA070F40E3585ULL, 0xD6990624D192E819ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG3, MSG2, 4);
    MSG0 = _mm256_add_epi32(MSG0, TMP);
    MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

    /* Rounds 48-51 */
    MSG = _mm256_add_epi32(MSG0, _mm256_set_epi64x(0x34B0BCB52748774CULL, 0x1E376C0819A4C116ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG0, MSG3, 4);
    MSG1 = _mm256_add_epi32(MSG1, TMP);
    MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
    MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

    /* Rounds 52-55 */
    MSG = _mm256_add_epi32(MSG1, _mm256_set_epi64x(0x682E6FF35B9CCA4FULL, 0x4ED8AA4A391C0CB3ULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG1, MSG0, 4);
    MSG2 = _mm256_add_epi32(MSG2, TMP);
    MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

    /* Rounds 56-59 */
    MSG = _mm256_add_epi32(MSG2, _mm256_set_epi64x(0x8CC7020884C87814ULL, 0x78A5636F748F82EEULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    TMP = _mm256_alignr_epi8(MSG2, MSG1, 4);
    MSG3 = _mm256_add_epi32(MSG3, TMP);
    MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

    /* Rounds 60-63 */
    MSG = _mm256_add_epi32(MSG3, _mm256_set_epi64x(0xC67178F2BEF9A3F7ULL, 0xA4506CEB90BEFFFAULL));
    STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
    MSG = _mm256_shuffle_epi32(MSG, 0x0E);
    STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

    /* Combine state  */
    STATE0 = _mm256_add_epi32(STATE0, ABEF_SAVE);
    STATE1 = _mm256_add_epi32(STATE1, CDGH_SAVE);

  }

  TMP = _mm256_shuffle_epi32(STATE0, 0x1B);       /* FEBA */
  STATE1 = _mm256_shuffle_epi32(STATE1, 0xB1);    /* DCHG */
  STATE0 = _mm256_blend_epi16(TMP, STATE1, 0xF0); /* DCBA */
  STATE1 = _mm256_alignr_epi8(STATE1, TMP, 8);    /* ABEF */

  /* Save state */
  _mm_storeu_si128((__m256i*) &state[0], STATE0);
  _mm_storeu_si128((__m256i*) &state[4], STATE1);
}



#if defined(TEST_MAIN)

#include <stdio.h>
#include <string.h>
int main(int argc, char* argv[])
{
    /* empty message with padding */
    uint8_t message[64];
    memset(message, 0x00, sizeof(message));
    message[0] = 0x80;

    /* initial state */
    uint32_t state[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    sha256_process_x86(state, message, sizeof(message));

    const uint8_t b1 = (uint8_t)(state[0] >> 24);
    const uint8_t b2 = (uint8_t)(state[0] >> 16);
    const uint8_t b3 = (uint8_t)(state[0] >>  8);
    const uint8_t b4 = (uint8_t)(state[0] >>  0);
    const uint8_t b5 = (uint8_t)(state[1] >> 24);
    const uint8_t b6 = (uint8_t)(state[1] >> 16);
    const uint8_t b7 = (uint8_t)(state[1] >>  8);
    const uint8_t b8 = (uint8_t)(state[1] >>  0);

    /* e3b0c44298fc1c14... */
    printf("SHA256 hash of empty message: ");
    printf("%02X%02X%02X%02X%02X%02X%02X%02X...\n",
        b1, b2, b3, b4, b5, b6, b7, b8);

    int success = ((b1 == 0xE3) && (b2 == 0xB0) && (b3 == 0xC4) && (b4 == 0x42) &&
                    (b5 == 0x98) && (b6 == 0xFC) && (b7 == 0x1C) && (b8 == 0x14));

    if (success)
        printf("Success!\n");
    else
        printf("Failure!\n");

    return (success != 0 ? 0 : 1);
}

#endif
