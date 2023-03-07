// ---------------------------------------------------------------------------+
//                               Assumptions                                  |
// 1- First word of hash is enough to hold the distinguish point and nserver  |
// e.g. with sha256 32bits is enough.                                         |
// - if this doesn't apply to your case edit to_which_server function         |
// see type of variable ones and ... in phase i functions                     |
// 2- The digest bits are the first n bits of the output of the hash function |
// ---------------------------------------------------------------------------+
//                          TABLE OF CONTENTS                                 |
//                                                                            |
// 1- Hash function generic setup                                             |
// 2- Long message attack main parameters                                     |
// 3- Configurations affected by MPI                                          |
// 4- SIMD configurations                                                     |
// -- a discover avx vector length                                            |
// -- set types for the vector                                                |
// ---------------------------------------------------------------------------+
///  @todo fix naming 
// annoying convention: sometimes we use size a meaning in bits
// and sometimes in byt


#ifndef LONG_MESSAGE_CONFIG
#define LONG_MESSAGE_CONFIG
#include <stdint.h>
#include <stdio.h>
#include "numbers_shorthands.h"
#include "confg_math_func.h"


/*  huge pages */
#include <sys/mman.h>
/* 2 MiB page  */
#define HPAGE_SIZE (1 << 21)

/* 1 GiB page  */
#define GPAGE_SIZE (1 << 30)




// ---------------------------------------------------------------------------+
//                   1- Hash function generic setup                           |
// ---------------------------------------------------------------------------+

/* change it only when you change the function */
#define WORD_SIZE 4 // bytes = 32 bit
#define WORD_SIZE_BITS (WORD_SIZE*8) // bytes = 32 bit
#define WORD_TYPE u32
#define NWORDS_DIGEST 8
#define NWORDS_STATE 8
#define NWORDS_INPUT 16
#define HASH_INPUT_SIZE (NWORDS_INPUT * WORD_SIZE) /* in bytes */
#define HASH_STATE_SIZE (NWORDS_STATE * WORD_SIZE) /* in bytes */ 
#define FILE_NAME_MAX_LENGTH 256






// -------------------------------------------------------------------------+
//                   2- Long message attack main parameters                 |
// -------------------------------------------------------------------------+
// Apologies: only n ≡ 0 mod 8 is allowed. This is not a feature.

// Let N := n / 8
/* bytes i.e n := 8*N bits */
#define N 12
/* record the whole state after each each interval has passed */
#define INTERVAL (1LL<<30)

 /* store 2^L elements in the dictionary  */

#define L 43



#define L_IN_BYTES CEILING(L, 8) /* How many bytes to accommedate L */

// wlog: we can stick to  power of 2, then dictionary might reject some
#define NHASHES (1LL<<L) // How many hashes we'll send to all dictionaries?

// nbits are zero, tphis will be defined according to the send latency
#define DIFFICULTY 2
/* How many words are needed to accommedate N */
#define N_NWORDS_CEIL CEILING(N, WORD_SIZE)
/* since there migh be an empty space on the right side of the word, we shift */
/* the mask that of the distinguished point test to the left */
#define SHIFT_AMOUNT (N_NWORDS_CEIL * 32 - N * 8)
#define DIST_PT_MASK_UNSHIFTED ((1 << DIFFICULTY) - 1)
#define DIST_PT_MASK (DIST_PT_MASK_UNSHIFTED << SHIFT_AMOUNT)


/* we are not going to hold more than 32 bits in a dict entry */
#define MAX_VAL_SIZE 32 /* in bits */
#define MAX_VAL_SIZE_BYTES 8 

/* how big is our counter */
#define CTR_TYPE u64


// sanity check
#if N > WORD_SIZE*NWORDS_DIGEST
  #error "you are trying to attack more bits than the hash function provides"
#endif


// -------------------------------------------------------------------------+



// -------------------------------------------------------------------------+
//                     3- Configurations effected by MPI                    |
// -------------------------------------------------------------------------+
// lines tagged with @python will be edited by a script                     |
// -------------------------------------------------------------------------+


/* edit manually */
#define NSERVERS 1
#define LOG2_NSERVERS BITS_TO_REPRESENT(NSERVERS)
#define DEFINED_BITS (LOG2_NSERVERS + DIFFICULTY) // @todo check
/* we might ignore few bits due to ceiling  */
#define DEFINED_BYTES CEILING(DEFINED_BITS, 8)




// How large is the dictionary in a server
/* #define NSLOTS_MY_NODE 8462608000LL */


// 7*4 GiB
#define TOTAL_RAM 20850444000LL
#define NRECEIVERS_PER_NODE 1
#define NSLOTS_MY_NODE (TOTAL_RAM / (VAL_SIZE_BYTES*NRECEIVERS_PER_NODE))




/* If each server has a different memory size then their candidates will have */
/* a different false positive probability since #discarded bits vary.         */
/* The following define scale the number according to the false positive.     */
/* The main obstacle in writing a formula as a macro is that I don't know a   */
/* way to get the ram size  macros. maybe a python script will overwrite the  */
/* definintion below.                                                         */




// PLEASE DO NOT EDIT// 
#ifndef LONG_MESSAGE_MPI_CONFIG

#define LONG_MESSAGE_MPI_CONFIG

#define TAG_DICT_SND 0 /* send digest tag */
/* rehashing defined part of long messageg is done for process i */
#define TAG_DONE_HASHING 1
#define TAG_RANDOM_MESSAGE 2
#define TAG_SND_DGST 3
#define TAG_MESSAGES_CANDIDATES 4
#define ARCHIVE_SERVER NSERVERS
// WHAT IS BUFF_SIZE? It doesn't seem to be used!
//#define BUFF_SIZE 1000  // holds `BUFF_SIZE` elements. @by_hand
#define PROCESS_QUOTA 1000000LL // i.e. send 10 digests to each server @by_hand

#endif // LONG_MESSAGE_MPI_CONFIG


/* #if LOG2_NSERVERS > N */
/*   #error "We have N < DEFINED_BYTES, this program is a bit stupid" */
/* #endif */













// -------------------------------------------------------------------------+
//                           3- SIMD configurations                         |
// -------------------------------------------------------------------------+
//                           Part a: Find Vector Size                       |
// -------------------------------------------------------------------------+
#ifdef __AVX512F__ /* F for Foundation */
#define REG_TYPE __m512i
#define REG_CMP_TYPE __mmask32
#define AVX_SIZE  512
#define ALIGNMENT 64
// Number of elements per simd register
#define SIMD_LEN  (AVX_SIZE / (8*VAL_SIZE_BYTES))




#define SIMD_LOAD_SI _mm512_load_si512
/* #define SIMD_TEST _mm512_test_epi64_mask */
#define SIMD_TEST is_not_zero
#define SIMD_SET1_EPI32 _mm512_set1_epi32
#define SIMD_SET1_EPI16 _mm512_set1_epi16
#define SIMD_SETZERO_SI  _mm512_setzero_si512

  
// todo check these 
#define SIMD_CMP_EPI32(X, Y) _mm512_cmpeq_epi32_mask(X, Y)
#define SIMD_CMP_EPI16(X, Y) _mm512_cmpeq_epi16_mask(X, Y)

#define SIMD_AND_EPI32 _mm512_and_epi32


#define SIMD_SET1_EPI8 _mm_set1_epi8
#define SIMD_CMP_EPI8 _mm_cmpeq_epi8

#elif defined(__AVX2__)
#define REG_TYPE __m256i
#define REG_CMP_TYPE __m256i
#define AVX_SIZE  256
#define ALIGNMENT 32
#define SIMD_LEN  (AVX_SIZE / (8*VAL_SIZE_BYTES))


#define SIMD_LOAD_SI _mm256_load_si256
#define SIMD_TEST _mm256_testz_si256


#define SIMD_SET1_EPI32 _mm256_set1_epi32
#define SIMD_SET1_EPI16 _mm256_set1_epi16
#define SIMD_SETZERO_SI _mm256_setzero_si256

#define SIMD_CMP_EPI32(X, Y) _mm256_movemask_epi8(_mm256_cmpeq_epi32(X, Y))
#define SIMD_CMP_EPI16(X, Y) _mm256_movemask_epi8(_mm256_cmpeq_epi16(X, Y))
#define SIMD_AND_EPI32  _mm256_andnot_si256 //@todo check!

// we better avoid this case :(
/* #define SIMD_SET1_EPI8 _mm256_set1_epi8 */
/* #define SIMD_CMP_EPI8 _mm256_cmpeq_epi8 */

/* #elif defined(__SSE4_1__) */
/* #define REG_TYPE __m128i */
/* #define REG_TYPE_CMP __m128i */
/* #define AVX_SIZE  128 */
/* #define ALIGNMENT 16 */
/* #define SIMD_LEN  (AVX_SIZE / (8*VAL_SIZE_BYTES)) */

/* #define SIMD_LOAD_SI _mm_load_si128 */
/* #define SIMD_TEST _mm_testz_si128 /\* as fast as RER C *\/ */


/* #define SIMD_SET1_EPI32 _mm_set1_epi32 */
/* #define SIMD_SET1_EPI16 _mm_set1_epi16 */

/* #define SIMD_CMP_EPI32 _mm_cmpeq_epi32 */
/* #define SIMD_CMP_EPI16 _mm_cmpeq_epi16 */

/* #define SIMD_SET1_EPI8 _mm_set1_epi8 */
/* #define SIMD_CMP_EPI8 _mm_cmpeq_epi8 */

#else
 #error "Please consider buying a modern device"
#endif // end avx testing

#define NHASH_LANES ((AVX_SIZE)/(WORD_SIZE_BITS))

// -------------------------------------------------------------------------+
//                     Part b: Set the best vector type                     |
// -------------------------------------------------------------------------+

// @todo check if the below comments need to be removed 

// Assumptions that need to be changed manually

// Dictionary Configurations:
// Recall: input -> store_as_idx, store_as_val
// we limit |store_as_val| <= 32bit

//  What is the minimum dictionary size we will use?
//#define L 32 /* i.e. it can hold a dictionary with nslots = 2^L */


// because the smallest step is one byte, we will not bother with shifting
// we will discard these bits. The program will work.
// e.g. imagine L = 34, then we have
// B0, B1, ..., B7, B8 mod 2^L will neglect 6 bits

// #define IDX_DISCARDED_BITS (L % 8) /* Values */
// Define VAL_SIZE_BYTES and decide what is the best type to use
/* Every thing will be stored as idx */
// #if N*8 <= L this case
//



// define intrinsics based on the the available max register size
// @todo start here 
#if N <= L_IN_BYTES
  #pragma message ("N is smaller than L!")
  #error "The code is not flexible to store every bit as index!"
#endif


// Having large valtype will make index smaller thus the filling rate
// will drop significantly

#if L_IN_BYTES <= N - DEFINED_BYTES - 4

#define VAL_SIZE_BYTES 4 /* byte */
#define VAL_TYPE u32 /* unsigned char */
// SIMD instructions has to be adapted according to the size
#define SIMD_SET1_VALTYPE SIMD_SET1_EPI32
#define SIMD_CMP_VALTYPE(X, Y) SIMD_CMP_EPI32(X, Y)

#else
#define VAL_SIZE_BYTES 2 /* byte */
#define VAL_TYPE u16 /* unsigned char */
// SIMD instructions has to be adapted according to the size
#define SIMD_SET1_VALTYPE SIMD_SET1_EPI16
#define SIMD_CMP_VALTYPE(X, Y) SIMD_CMP_EPI16(X, Y)

#endif



/* [distinguished point bits || nserver || y || idx || val || z ] */
// The discarded bits are:
// 1- y since (distinguished points || nserver ) may not be multiple of 8bits
// 2- z

// @todo this is not satisfactory!
#define DISCARDED_BITS MAX((8 * N - L - 8 * VAL_SIZE_BYTES), 0)
/* #define DISCARDED_BITS MAX((8 * N - L - 8 * VAL_SIZE_BYTES - DEFINED_BITS), 0) */

//#define DISCARDED_BITS (8*N - L - 8 * VAL_SIZE_BYTES)
// We need 2^#disacrded_bits candidates, we expect each server generate

#define NNEEDED_CND MAX(( (1LL << (DISCARDED_BITS+3)) ),\
			  1)   /* @python  */

/* #define NNEEDED_CND 5000 */
		

// (2^#disacrded_bits) / NSERVERS, however, it's not a strict requirement.
//#define NNEEDED_CND_THIS_SERVER MAX(((1LL << DISCARDED_BITS) >> NSERVERS), 1)
  
// a candidate is a hash that was found in the dictionary.
// since we may have false positive, we need to get generat number of them





#endif




