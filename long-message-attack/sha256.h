/*********************************************************************
* Filename:   sha256.h
* Author:     Brad Conte (brad AT bradconte.com)
* Copyright:
* Disclaimer: This code is presented "as is" without any guarantees.
* Details:    Defines the API for the corresponding SHA1 implementation.
* Notes: this code diverged from the original code 
*********************************************************************/

#ifndef SHA256_H
#define SHA256_H

/*************************** HEADER FILES ***************************/
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "types.h" // digest type
#include "dict.h" // dict type 
/****************************** MACROS ******************************/
#define SHA256_BLOCK_SIZE 32            // SHA256 outputs a 32 byte digest
// for the modified version @ahmed
#define OUTPUT_SIZE_BITS 256 // this has to be a parameter not a macro :)

/**************************** DATA TYPES ****************************/
typedef unsigned char BYTE;             // 8-bit byte
typedef unsigned int  WORD;             // 32-bit word, change to "long" for 16-bit machines
//@ahmed since the output will be tunrcate we migh opt for a small type
typedef unsigned int OUTPUT_TYPE; // control the output type

typedef struct {
	BYTE data[64];
	WORD datalen;
	unsigned long long bitlen;
	WORD state[8];
} SHA256_CTX;


// todo remove this and use digest union 
typedef union {
  WORD* state;
  char* sate_as_bytes;
} STATE;

/*********************** FUNCTION DECLARATIONS **********************/
void sha256_init(SHA256_CTX *ctx);
/* void sha256_update(SHA256_CTX *ctx, const BYTE data[], BYTE** intermediate, size_t len, int output_size_bits); */
/* void sha256_final(SHA256_CTX *ctx, BYTE hash[], int output_size_bits); */
void sha256_transform(SHA256_CTX *ctx, const BYTE data[], int output_size_bits);
void sha256_update(SHA256_CTX *ctx, const BYTE data[], size_t len,
		   int output_size_bits, dict* d, void (*add_to_dict)(dict*, char*, size_t, size_t));
void sha256_final(SHA256_CTX *ctx, BYTE hash[], int output_size_bits);
void print_intermediate(SHA256_CTX* ctx);
#endif   // SHA256_H
