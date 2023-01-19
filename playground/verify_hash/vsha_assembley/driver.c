#include <stdio.h>
#include <inttypes.h>

typedef uint32_t u32;
typedef uint8_t u8;

#include "arch_avx2_type1.h"
#define SHA256_H0 0x6a09e667
#define SHA256_H1 0xbb67ae85
#define SHA256_H2 0x3c6ef372
#define SHA256_H3 0xa54ff53a
#define SHA256_H4 0x510e527f
#define SHA256_H5 0x9b05688c
#define SHA256_H6 0x1f83d9ab
#define SHA256_H7 0x5be0cd19

void sha256_mb_init_digest(uint32_t *digest)
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

// sha_mb_generic_init(state->args.digest

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