/// we define here how large uint

#ifndef LONG_MESSAGE_CONFIG
#define LONG_MESSAGE_CONFIG
#include <stdint.h>
#include "numbers_shorthands.h"


// hash function specifications
// Note: we have not defined that attack parametes
// since we use sha256
#define WORD_SIZE 4 // bytes = 32 bit  
#define NWORDS_DIGEST 8 
#define NWORDS_INPUT 16


// -------------------------------------------------------------------------+
//                     Long message attack main parameters                  |
// -------------------------------------------------------------------------+
// Apologies: only n ≡ 0 mod 8 is allowed. This is not a feature.


// Let N := n / 8.  
#define N 12 // bytes i.e n = 8*N bits
#define DIFFICULTY 4 // bits are zero

// Assumptions that need to be changed manually

// Dictionary Configurations:
// Recall: input -> store_as_idx, store_as_val
// we limit |store_as_val| <= 32bit

//  What is the minimum dictionary size we will use?
#define L_MIN 32 /* i.e. it can hold a dictionary with nslots = 2^L_MIN */
#define MAX_VAL_SIZE 32 // we are not going to hold more than 32 bits

// because the smallest step is one byte, we will not bother with shifting
// we will discard these bits. The program will work.
// e.g. imagine L_MIN = 34, then we have
// B0, B1, ..., B7, B8 mod 2^L_MIN will neglect 6 bits

#define IDX_DISCARDED_BITS (L_MIN % 8) /* Values */
// Define VAL_SIZE and decide what is the best type to use
/* Every thing will be stored as idx */
// #if N*8 <= L_MIN this case

#if N * 8 <= L_MIN
  #pragma message ("Please decrease MAX_VAL_SIZE")
  #error "The code is not flexible to store every bit as index!"

// we store everything as index except one byte
//+ @todo when changing MAX_VAL_SIZE adopt the conditions below
#elif (N * 8 < L_MIN + 8) && (N*8 >= L_MIN) 
  #define VAL_SIZE 1 /* byte */
  #define VAL_TYPE u8  /* unsigned char */
  /* we have to get more positive probes depending on the next value  */
  #define DISCARDED_BITS IDX_DISCARDED_BITS + 8*(N - VAL_SIZE)

#elif (N * 8 < L_MIN + 16) && (N*8 >= L_MIN + 8) 
  #define VAL_SIZE 2 /* byte */
  #define VAL_TYPE u16 /* unsigned char */

// Typical condition when in terms of MAX_VAL_SIZE
#elif (N * 8 < L_MIN + MAX_VAL_SIZE) && (N*8 >= L_MIN + (MAX_VAL_SIZE/2) )
  #define VAL_SIZE 4 /* byte */
  #define VAL_TYPE u32 /* unsigned char */


#else /* The MAX_VAL_SIZE */
  #define VAL_SIZE 4 /* byte */
  #define VAL_TYPE u32 /* unsigned char */

#endif // define VAL_SIZE
// -------------------------------------------------------------------------+

 
// depending on avx register length, on my laptop 256
//  Change them together, @todo write a code to automate writing these values
// -------------------------------------------- +
#define AVX_SIZE  256                        // |
#define ALIGNMENT 32                         // |
// Number of elements per simd register         |
#define SIMD_LEN  (AVX_SIZE / (8*VAL_SIZE))  // |
// -------------------------------------------- +




#ifndef LONG_MESSAGE_MPI_CONFIG
#define LONG_MESSAGE_MPI_CONFIG
// MPI configurations
// #define NSERVERS 10 we can get it from mpi_comm_size
#define LOG2_NSERVERS 4 // = ceil(log2(NSERVERS))
#define BUFF_SIZE 1000  // holds `BUFF_SIZE` elements.

#define PROCESS_QUOTA 100 // i.e. send 10 digests to each server

#define NWORDS_OFFSET 4 // use 128 bits as offsets to find message
#endif // LONG_MESSAGE_MPI_CONFIG







#endif




