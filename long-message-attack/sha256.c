/*********************************************************************
* Filename:   sha256.c
* Author:     Brad Conte (brad AT bradconte.com)
* Copyright:
* Disclaimer: This code is presented "as is" without any guarantees.
* Details:    Implementation of the SHA-256 hashing algorithm.
              SHA-256 is one of the three algorithms in the SHA2
              specification. The others, SHA-384 and SHA-512, are not
              offered in this implementation.
              Algorithm specification can be found here:
               * http://csrc.nist.gov/publications/fips/fips180-2/fips180-2withchangenotice.pdf
              This implementation uses little endian byte order.
*
* Changes by @ahmed:
* - removed merkle-damgard strngthen
* - truncate the digest
* - intermediate states are stored in a dictionary
* - states are represented as an instance of union
*   not sure how the complex the code is becoming :(
*********************************************************************/

/*************************** HEADER FILES ***************************/

#include <stdint.h>
#include <stdlib.h>
#include <memory.h>
#include <sys/types.h>
#include "dict.h"
#include "sha256.h"
#include "types.h"

/****************************** MACROS ******************************/
#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

/**************************** VARIABLES *****************************/
static const WORD k[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

/*********************** FUNCTION DEFINITIONS ***********************/
void print_intermediate(SHA256_CTX* ctx){
  for (int i=0; i<8; i++)
    printf("0x%08x, ", ctx->state[i]);
  puts("");
}

void truncate(SHA256_CTX* ctx, int output_size_bits){
  // ctx->state[i] is a 32 bit block
  // we order the bits of state[0], state[1], ..., etc as:
  // b0 b1 ... b31    b32 b33 ... b63   b64 b65 ... b95
  //   state[0]           state[1]       state[3] .....
  // Goal: set all bits bi s.t. i>= output_size_bits to 0
  // e.g. output_size_bits = 4
  // b0 b1 b2 b3 0 ...0

  // all blocks after ceil(output_size_bits/32) should be zero
  // For the singular case when someone choses the output to be zero
  assert(output_size_bits != 0);
  // printf("output_size=%d\n", output_size_bits);
  // see defintion of `rem` below. the extra term is to ensure that ceil adds one
  // the division is exact which mean the last active block has all bits active
  int i =  ceil((float)output_size_bits/32  + ((double) 2 / 32));
  // printf("we are going to truncate from i=%d\n", i);

  for (int j = i; j<8; ++j){
    // printf("ctx->state[%d] was %8x \n", j, ctx->state[i]);
    ctx->state[j] = 0;
    // printf("ctx->state[%d] is %8x \n", j, ctx->state[i]);
  }
  // and with output_size_bits ones
  // this convluted formula to deal with the case when all the bits of the last should be active
  // the mod 32 doesn't caputre 32 bits number
  // another solution to add ε < (1/32) inside the ceil within the definiton of i
  uint64_t rem = output_size_bits % 32;
  // printf("rem=%lu bits\n", rem);
  rem = ((uint64_t) 1<<rem) - 1; // 111...1 /
  ctx->state[i - 1] = ctx->state[i-1] & rem;

}



void sha256_transform(SHA256_CTX *ctx, const BYTE data[], int output_size_bits)
{
	WORD a, b, c, d, e, f, g, h, i, j, t1, t2, m[64];

	for (i = 0, j = 0; i < 16; ++i, j += 4)
		m[i] = (data[j] << 24) | (data[j + 1] << 16) | (data[j + 2] << 8) | (data[j + 3]);
	for ( ; i < 64; ++i)
		m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

	a = ctx->state[0];
	b = ctx->state[1];
	c = ctx->state[2];
	d = ctx->state[3];
	e = ctx->state[4];
	f = ctx->state[5];
	g = ctx->state[6];
	h = ctx->state[7];

	for (i = 0; i < 64; ++i) {
		t1 = h + EP1(e) + CH(e,f,g) + k[i] + m[i];
		t2 = EP0(a) + MAJ(a,b,c);
		h = g;
		g = f;
		f = e;
		e = d + t1;
		d = c;
		c = b;
		b = a;
		a = t1 + t2;
	}

	ctx->state[0] += a;
	ctx->state[1] += b;
	ctx->state[2] += c;
	ctx->state[3] += d;
	ctx->state[4] += e;
	ctx->state[5] += f;
	ctx->state[6] += g;
	ctx->state[7] += h;

	/* puts("before truncation"); */
	/* print_intermediate(ctx); */
	// @note we are going truncate the internal
	// state rather we will truncate outside the
	// sha256 transform
	// truncate(ctx, output_size_bits);
}

void sha256_init(SHA256_CTX *ctx)
{
        
	ctx->datalen = 0;
	ctx->bitlen = 0;
	ctx->state[0] = 0x6a09e667;
	ctx->state[1] = 0xbb67ae85;
	ctx->state[2] = 0x3c6ef372;
	ctx->state[3] = 0xa54ff53a;
	ctx->state[4] = 0x510e527f;
	ctx->state[5] = 0x9b05688c;
	ctx->state[6] = 0x1f83d9ab;
	ctx->state[7] = 0x5be0cd19;

	
}


// void sha256_update(SHA256_CTX *ctx, const BYTE data[], BYTE** intermediate,
// size_t len, int output_size_bits)
void sha256_update(SHA256_CTX *ctx, const BYTE data[], size_t len, int output_size_bits)//, dict* d, void (*add_to_dict)(dict*, dict_key*, size_t, size_t))

{
        
        // printf("update len = %lu bytes\n", len);
        WORD i;	
	//  ctx->state;
	// int output_size_bytes = (int) ceil((double) output_size_bits/8);
	// printf("output size in bytes=%d\n", output_size_bytes);
	for (i = 0; i < len; ++i) {
		ctx->data[ctx->datalen] = data[i];
		ctx->datalen++;

		if (ctx->datalen == 64) {
		  sha256_transform(ctx, ctx->data, output_size_bits);
		  ctx->bitlen += 512;
		  ctx->datalen = 0;
		  // @ahmed check
		  
		  // printf("after %lld bits\n", ctx->bitlen);
		  // print_intermediate(ctx);
		  // add_element_to_dictionary(dict *dictionary, char *key, size_t value, size_t input_size)

		  //		  add_to_dict(d, (dict_key *) ctx->state,((int) i/64 )+ 1, output_size_bytes );
		  if (is_there_duplicate)
		    break;
		  
		  //puts("");
		}
	}
}

// void sha256_final(SHA256_CTX *ctx, BYTE hash[], int output_size_bits)
void sha256_final(SHA256_CTX *ctx, BYTE hash[], int output_size_bits)
{
        // WORD i;

	// i = ctx->datalen;
	// Pad whatever data is left in the buffer.
	// this padding has to be edited
	// ctx -> datalen is in number of bytes
	// 10* pad @ahmed
	if (ctx->datalen < 64){
	  ctx->data[ctx->datalen++] = 0x80;
	}
	while (ctx->datalen < 64){
	  ctx->data[ctx->datalen++] = 0x00;
	}



	/* if (ctx->datalen < 56) { */
	/* 	ctx->data[i++] = 0x80; */
	/* 	while (i < 56) */
	/* 		ctx->data[i++] = 0x00; */
	/* } */
	/* else { */
	/* 	ctx->data[i++] = 0x80; */
	/* 	while (i < 64) */
	/* 		ctx->data[i++] = 0x00; */
	/* 	sha256_transform(ctx, ctx->data); */
	/* 	memset(ctx->data, 0, 56); */
	/* } */

	// Append to the padding the total message's length in bits and transform.
	ctx->bitlen += ctx->datalen * 8;
	// this the strengthen of Merkle-Damgard, ignore it @ahmed
        /* ctx->data[63] = ctx->bitlen; */
	/* ctx->data[62] = ctx->bitlen >> 8; */
	/* ctx->data[61] = ctx->bitlen >> 16; */
	/* ctx->data[60] = ctx->bitlen >> 24; */
	/* ctx->data[59] = ctx->bitlen >> 32; */
	/* ctx->data[58] = ctx->bitlen >> 40; */
	/* ctx->data[57] = ctx->bitlen >> 48; */
	/* ctx->data[56] = ctx->bitlen >> 56; */
	sha256_transform(ctx, ctx->data, output_size_bits);

	// for now let's forget the big-endian conversion @ahmed
	/* // Since this implementation uses little endian byte ordering and SHA uses big endian, */
	/* // reverse all the bytes when copying the final state to the output hash. */
	/* for (WORD i = 0; i < 4; ++i) { */
	/* 	hash[i]      = (ctx->state[0] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 4]  = (ctx->state[1] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 8]  = (ctx->state[2] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 20] = (ctx->state[5] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 24] = (ctx->state[6] >> (24 - i * 8)) & 0x000000ff; */
	/* 	hash[i + 28] = (ctx->state[7] >> (24 - i * 8)) & 0x000000ff; */
	/* } */
}


// @ahmed
void truncate_state_get_digest(uint64_t* dst, SHA256_CTX* ctx, int n_of_bits){
  /// We extract the digest from ctx and save it in dst
  // uint64_t dst[2] is fixed now // 128 bits, this is a limitation
  // it should be 256 for sha256 :)
  // n_of_bits is how many bits is the compression function output
  

  dst[0] = ctx->state[0] + (((uint64_t) ctx->state[1])<<32);
  if (n_of_bits < 64){
    uint64_t ones =   (((uint64_t) 1) << (n_of_bits)) - 1;
    dst[0] = dst[0] & ones;
    dst[1] = 0;

    #ifdef VERBOSE_LEVEL
    printf("ones=%lu\n", ones);
    printf("state[0]=%x, state[1]=%x\n", ctx->state[0], ctx->state[1]);
    puts("");
    #endif // VERBOSE_LEVEL

    
  } else if (n_of_bits < 128) {
    // since the number of bits is higher or equal 64
    // we need to work on the second element of the dst
    n_of_bits = n_of_bits - 64;

    #ifdef VERBOSE_LEVEL
    uint64_t ones =  ( ((uint64_t) 1<<(n_of_bits)) - 1);
    printf("ones=%lu, n_of_bits=%d\n", ones, n_of_bits);
    printf("state[0]=%x, state[1]=%x\n", ctx->state[0], ctx->state[1]);
    printf("state[2]=%x, state[2]=%x\n", ctx->state[2], ctx->state[3]);
    puts("");
    #endif // VERBOSE_LEVEL

    // copy 64bits from state
    dst[1] = ctx->state[2] + (((uint64_t) ctx->state[3])<<32);
    // truncate it if necessary
    dst[1] = dst[1] & ( ((uint64_t) 1<<(n_of_bits)) - 1);
  } else { // 128 bits limit
    dst[1] = ctx->state[2] + (((uint64_t) ctx->state[3])<<32);
  }
    
  
}
