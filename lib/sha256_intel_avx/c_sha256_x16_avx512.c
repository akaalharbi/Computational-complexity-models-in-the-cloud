#include <stdio.h>
#include <stdint.h>
#include "c_sha256_x16_avx512.h"
#define AVX2_NLANES_SHA256 16

typedef uint32_t u32;
typedef uint8_t u8;

#include "arch_avx512_type1.h"
#define SHA256_H0 0x6a09e667
#define SHA256_H1 0xbb67ae85
#define SHA256_H2 0x3c6ef372
#define SHA256_H3 0xa54ff53a
#define SHA256_H4 0x510e527f
#define SHA256_H5 0x9b05688c
#define SHA256_H6 0x1f83d9ab
#define SHA256_H7 0x5be0cd19




static void sha256_mb_init_digest_avx512(uint32_t *digest)
{
	/* 8 is sufficient for AVX2, 16 goes all the way to AVX512 */
	for (int lane = 0; lane < 16; lane++) {
        digest[lane + 0*16] = SHA256_H0;
        digest[lane + 1*16] = SHA256_H1;
        digest[lane + 2*16] = SHA256_H2;
        digest[lane + 3*16] = SHA256_H3;
        digest[lane + 4*16] = SHA256_H4;
        digest[lane + 5*16] = SHA256_H5;
        digest[lane + 6*16] = SHA256_H6;
        digest[lane + 7*16] = SHA256_H7;
    }
}


static void sha256_init_digest_avx512(uint32_t *digest, uint32_t tr_states[16*8])
{
	/* 8 is sufficient for AVX2, 16 goes all the way to AVX512 */
	for (int i = 0; i < 16*8; ++i) {
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
        digest[i] = tr_states[i];
    }
}


uint32_t* sha256_multiple_x16(uint8_t msg[16][64]){
  //---------------------------------------------------------------------------+
  // The unfortunate convention 16 is the number of lanes, in avx2 it will only|
  // use 8 of them. However, it's mandatory to supply 16 messages.             |
  // NOTE: the returned digest is transposed it's important to transpose it    |
  //       to get a correct results.                                           |
  // e.g. ith message digest can be found by reading:                          |
  // args.digest[i + 16*0]                                                     |
  // args.digest[i + 16*1]                                                     |
  // ...                                                                       |
  // args.digest[i + 16*7]                                                     |  
  //---------------------------------------------------------------------------+
  
  // 32 bytes = 512 bits (input size)
  static SHA256_ARGS args; /* test static */
  sha256_mb_init_digest_avx512(args.digest);

  for (int lane=0; lane<AVX2_NLANES_SHA256; ++lane) {
    args.data_ptr[lane] = msg[lane];
  }
  call_sha256_x16_avx512_from_c(&args, 1);


  return args.digest;
};

uint32_t *sha256_multiple_x16_tr(uint8_t msg[16][64], uint32_t tr_states[16 * 8])
{

  /* this function takes data NON-TRASNPOSED and states TRANSPOSED */
  
  // 32 bytes = 512 bits (input size)
  static SHA256_ARGS args; /* test static */
  /* if the function was called before don't init the data*/
  static int inited = 0;

  if (!inited) {
    sha256_init_digest_avx512(args.digest, tr_states);
  }


  for (int lane=0; lane<AVX2_NLANES_SHA256; ++lane) {
    args.data_ptr[lane] = msg[lane];
  }
  call_sha256_x16_avx512_from_c(&args, 1);
  inited = 1;

  return args.digest;
  
}  




#ifdef TESTMAIN
int main()
{

	/* test with msg = 000.........0000 */
	{
		SHA256_ARGS args;
		sha256_mb_init_digest(args.digest);
		
		/* 512 bit = 16 x u32 */
		u32 msg[16] = {};
		for (int lane = 0; lane < 16; lane++)
			args.data_ptr[lane] = (u8 *) msg;
		call_sha256_oct_avx2_from_c(&args, 1);
		
		for (int i = 0; i < 16; i++) {
			printf("lane %d: ", i);
			for (int j = 0; j < 8; j++)
				printf("%08x ", args.digest[i + 16*j]);
			printf("\n");
		}
	}

	/* test with msg = 00 01 02 ......... 3e 3f + per-lane modification */
	{
		printf("===============================================\n");
		SHA256_ARGS args;
		sha256_mb_init_digest(args.digest);
		/* 512 bit = 16 x u32 */
		u32 msg[16][16] = {};
		for (int lane = 0; lane < 16; lane++)
			args.data_ptr[lane] = (u8 *) msg[lane];
				
		for (int i = 0; i < 8; i++) {
			msg[i][0] = 0x10101010 * i;
			msg[i][1] = 0x10101010 * i;
			msg[i][2] = 0x10101010 * i;
			msg[i][3] = 0x10101010 * i;
			msg[i][4] = 0x10101010 * i;
			msg[i][5] = 0x10101010 * i;
			msg[i][6] = 0x10101010 * i;
			msg[i][7] = 0x10101010 * i;
			msg[i][8] = 0x10101010 * i;
			msg[i][9] = 0x10101010 * i;
			msg[i][10] = 0x10101010 * i;
			msg[i][11] = 0x10101010 * i;
			msg[i][12] = 0x10101010 * i;
			msg[i][13] = 0x10101010 * i;
			msg[i][14] = 0x10101010 * i;
			msg[i][15] = 0x10101010 * i;
		}	
		call_sha256_oct_avx2_from_c(&args, 1);
		for (int i = 0; i < 16; i++) {
			printf("lane %d: ", i);
			for (int j = 0; j < 8; j++)
				printf("%08x ", args.digest[i + 16*j]);
			printf("\n");
		}
	}
}

#endif
