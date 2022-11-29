/// we define here how large uint

#ifndef LONG_MESSAGE_CONFIG
#define LONG_MESSAGE_CONFIG
#include <stdint.h>

// dead weight they have to be reomved
#define UINT uint64_t // this
#define BYTE unsigned char



// hash function specific
// since we use sha256
#define WORD_SIZE 4 // bytes = 32 bit  
#define NWORDS_DIGEST 8 
#define NWORDS_INPUT 16 


// depending on avx register length, on my laptop 256
//  Change them together, @todo write a code to automate writing these values
// ---------------------- |
#define AVX_SIZE  256 //  |
#define ALIGNMENT 32  //  |
// Number of elements per simd register
#define SIMD_LEN  8   //  |
// ---------------------- |


// Long message attack tuning:
// Difficulty level
#define DIFFICULTY 4 // bits are zero


#ifndef LONG_MESSAGE_MPI_CONFIG
#define LONG_MESSAGE_MPI_CONFIG
// MPI configurations
#define NSERVERS 10 
#define LOG2_NSERVERS 4 // = ceil(log2(NSERVERS))
#define BUFF_SIZE 100  // holds `BUFF_SIZE` elements.

#define MY_QUOTA 10 // i.e. send 10 digests to each server

#define NWORDS_OFFSET 4 // use 128 bits as offsets to find message
#endif // LONG_MESSAGE_MPI_CONFIG







#endif




