
// #include <inttypes.h>
#include <stdio.h>
#include <stdint.h>
#include "vsha256.h"


/*
 * /!\ this code does not care about endianness. Some conversions might be in order.
 */



/*
 * This implementation of SHA256 is taken almost straight out of openssl,
 * with minor modifications.  It is OBLIVIOUS of vector length. This relies
 * on gcc's builtin support for vectors. This code is not portable to other
 * compilers.
 *
 * The same code would work with AVX-512, with only minor modifications.
 * By Charles Bouilaguet (I hope the spelling is correct)
 */


/*
 * vectorized sha256 implementation. This operates on the "vector" type, whatever this is.
 */

vector vbroadcast(u32 x)
{
	vector b = { x, x, x, x, x, x, x, x };
	return b;
}

void __sha256_init(vector state[8])
{
	state[0] = vbroadcast(0x6a09e667);
	state[1] = vbroadcast(0xbb67ae85);
	state[2] = vbroadcast(0x3c6ef372);
	state[3] = vbroadcast(0xa54ff53a);
	state[4] = vbroadcast(0x510e527f);
	state[5] = vbroadcast(0x9b05688c);
	state[6] = vbroadcast(0x1f83d9ab);
	state[7] = vbroadcast(0x5be0cd19);
}

static const u32 K256[64] = {
	0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
	0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
	0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
	0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
	0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
	0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
	0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
	0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
	0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
	0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
	0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
	0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
	0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
	0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
	0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
	0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

vector vK256[64];

void vsha256_setup()
{
	for (int i = 0; i < 64; i++)
		vK256[i] = vbroadcast(K256[i]);
}

/*
 * FIPS specification refers to right rotations, while our ROTATE macro
 * is left one. This is why you might notice that rotation coefficients
 * differ from those observed in FIPS document by 32-N...
 */
vector static inline ROTATE(vector a, int n)
{
	return (((a) << (n)) | (((a) & 0xffffffff) >> (32 - (n))));
}

vector static inline Sigma0(vector x)
{
	return (ROTATE((x), 30) ^ ROTATE((x), 19) ^ ROTATE((x), 10));
}

vector static inline Sigma1(vector x)
{
	return (ROTATE((x), 26) ^ ROTATE((x), 21) ^ ROTATE((x), 7));
}

vector static inline sigma0(vector x)
{
	return (ROTATE((x), 25) ^ ROTATE((x), 14) ^ ((x) >> 3));
}

vector static inline sigma1(vector x)
{
	return (ROTATE((x), 15) ^ ROTATE((x), 13) ^ ((x) >> 10));
}

vector static inline Ch(vector x, vector y, vector z)
{
	return (((x) & (y)) ^ ((~(x)) & (z)));
}

vector static inline Maj(vector x, vector y, vector z)
{
	return (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)));
}

void __sha256_transform(vector state[], const vector W[])
{
	vector X[16];
	vector a = state[0];
	vector b = state[1];
	vector c = state[2];
	vector d = state[3];
	vector e = state[4];
	vector f = state[5];
	vector g = state[6];
	vector h = state[7];

	/* first 16 rounds (consume the input message W as-is, copy it into X) */
	for (int i = 0; i < 16; i += 8) {
		X[0 + i] = W[0 + i];
		vector T10 = X[0 + i] + h + Sigma1(e) + Ch(e, f, g) + K256[0 + i];
		h = Sigma0(a) + Maj(a, b, c);
		d += T10;
		h += T10;

		X[1 + i] = W[1 + i];
		vector T11 = X[1 + i] + g + Sigma1(d) + Ch(d, e, f) + K256[1 + i];
		g = Sigma0(h) + Maj(h, a, b);
		c += T11;
		g += T11;

		X[2 + i] = W[2 + i];
		vector T12 = X[2 + i] + f + Sigma1(c) + Ch(c, d, e) + K256[2 + i];
		f = Sigma0(g) + Maj(g, h, a);
		b += T12;
		f += T12;

		X[3 + i] = W[3 + i];
		vector T13 = X[3 + i] + e + Sigma1(b) + Ch(b, c, d) + K256[3 + i];
		e = Sigma0(f) + Maj(f, g, h);
		a += T13;
		e += T13;

		X[4 + i] = W[4 + i];
		vector T14 = X[4 + i] + d + Sigma1(a) + Ch(a, b, c) + K256[4 + i];
		d = Sigma0(e) + Maj(e, f, g);
		h += T14;
		d += T14;

		X[5 + i] = W[5 + i];
		vector T15 = X[5 + i] + c + Sigma1(h) + Ch(h, a, b) + K256[5 + i];
		c = Sigma0(d) + Maj(d, e, f);
		g += T15;
		c += T15;

		X[6 + i] = W[6 + i];
		vector T16 = X[6 + i] + b + Sigma1(g) + Ch(g, h, a) + K256[6 + i];
		b = Sigma0(c) + Maj(c, d, e);
		f += T16;
		b += T16;

		X[7 + i] = W[7 + i];
		vector T17 = X[7 + i] + a + Sigma1(f) + Ch(f, g, h) + K256[7 + i];
		a = Sigma0(b) + Maj(b, c, d);
		e += T17;
		a += T17;
	}

	/* next 48 rounds (perform message expansion on the fly inside X) */
	for (int i = 16; i < 64; i += 8) {
		vector s00 = sigma0(X[(i + 1) & 0x0f]);
		vector s10 = sigma1(X[(i + 14) & 0x0f]);
		X[(i + 0) & 0x0f] += s00 + s10 + X[(i + 9) & 0x0f];
		vector T10 = X[(i + 0) & 0x0f] + h + Sigma1(e) + Ch(e, f, g) + K256[i + 0];
		h = Sigma0(a) + Maj(a, b, c);
		d += T10;
		h += T10;

		vector s01 = sigma0(X[(i + 2) & 0x0f]);
		vector s11 = sigma1(X[(i + 15) & 0x0f]);
		X[(i + 1) & 0x0f] += s01 + s11 + X[(i + 10) & 0x0f];
		vector T11 = X[(i + 1) & 0x0f] + g + Sigma1(d) + Ch(d, e, f) + K256[i + 1];
		g = Sigma0(h) + Maj(h, a, b);
		c += T11;
		g += T11;

		vector s02 = sigma0(X[(i + 3) & 0x0f]);
		vector s12 = sigma1(X[(i + 16) & 0x0f]);
		X[(i + 2) & 0x0f] += s02 + s12 + X[(i + 11) & 0x0f];
		vector T12 = X[(i + 2) & 0x0f] + f + Sigma1(c) + Ch(c, d, e) + K256[i + 2];
		f = Sigma0(g) + Maj(g, h, a);
		b += T12;
		f += T12;

		vector s03 = sigma0(X[(i + 4) & 0x0f]);
		vector s13 = sigma1(X[(i + 17) & 0x0f]);
		X[(i + 3) & 0x0f] += s03 + s13 + X[(i + 12) & 0x0f];
		vector T13 = X[(i + 3) & 0x0f] + e + Sigma1(b) + Ch(b, c, d) + K256[i + 3];
		e = Sigma0(f) + Maj(f, g, h);
		a += T13;
		e += T13;

		vector s04 = sigma0(X[(i + 5) & 0x0f]);
		vector s14 = sigma1(X[(i + 18) & 0x0f]);
		X[(i + 4) & 0x0f] += s04 + s14 + X[(i + 13) & 0x0f];
		vector T14 = X[(i + 4) & 0x0f];
		T14 += d + Sigma1(a) + Ch(a, b, c) + K256[i + 4];
		d = Sigma0(e) + Maj(e, f, g);
		h += T14;
		d += T14;

		vector s05 = sigma0(X[(i + 6) & 0x0f]);
		vector s15 = sigma1(X[(i + 19) & 0x0f]);
		X[(i + 5) & 0x0f] += s05 + s15 + X[(i + 14) & 0x0f];
		vector T15 = X[(i + 5) & 0x0f] + c + Sigma1(h) + Ch(h, a, b) + K256[i + 5];
		c = Sigma0(d) + Maj(d, e, f);
		g += T15;
		c += T15;

		vector s06 = sigma0(X[(i + 7) & 0x0f]);
		vector s16 = sigma1(X[(i + 20) & 0x0f]);
		X[(i + 6) & 0x0f] += s06 + s16 + X[(i + 15) & 0x0f];
		vector T16 = X[(i + 6) & 0x0f] + b + Sigma1(g) + Ch(g, h, a) + K256[i + 6];
		b = Sigma0(c) + Maj(c, d, e);
		f += T16;
		b += T16;

		vector s07 = sigma0(X[(i + 8) & 0x0f]);
		vector s17 = sigma1(X[(i + 21) & 0x0f]);
		X[(i + 7) & 0x0f] += s07 + s17 + X[(i + 16) & 0x0f];
		vector T17 = X[(i + 7) & 0x0f] + a + Sigma1(f) + Ch(f, g, h) + K256[i + 7];
		a = Sigma0(b) + Maj(b, c, d);
		e += T17;
		a += T17;
	}
	state[0] += a;
	state[1] += b;
	state[2] += c;
	state[3] += d;
	state[4] += e;
	state[5] += f;
	state[6] += g;
	state[7] += h;
}

/*****************************************************************************/

/*
 * Public interface. This encapsulates the "vector" type and exposes nicer functions
 * to the "end user" of the vectorized code.  This says "1 vector == 8 x u32".
 */
void vsha256_init(u32 state[8][8])
{
	__sha256_init((vector *) state);
}

void vsha256_transform(u32 state[8][8], const u32 W[16][8])
{
	__sha256_transform((vector *) state, (const vector *)W);
}

/* <----------- no "vector" beyond this point */

/*****************************************************************************/

/*
 * Example use
 */


#ifdef MAIN_VSHA256
int main()
{
	/* setup constants */
	vsha256_setup();

	/* test with msg = 000.........0000 */
	{
		u32 AVX_ALIGNED h[8][8];
		vsha256_init(h);
		u32 AVX_ALIGNED msg[16][8] = { };
		vsha256_transform(h, msg);
		for (int i = 0; i < 8; i++) {
			printf("lane %d: ", i);
			for (int j = 0; j < 8; j++)
				printf("%08x ", h[j][i]);
			printf("\n");
		}
	}

	/* test with msg = 00 01 02 ......... 3e 3f + per-lane modification */
	{
		printf("===============================================\n");
		u32 AVX_ALIGNED h[8][8];
		vsha256_init(h);
		u32 AVX_ALIGNED msg[16][8];
		for (int i = 0; i < 8; i++)
			msg[0][i] = 0x00010203 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[1][i] = 0x04050607 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[2][i] = 0x08090a0b + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[3][i] = 0x0c0d0e0f + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[4][i] = 0x10111213 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[5][i] = 0x14151617 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[6][i] = 0x18191a1b + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[7][i] = 0x1c1d1e1f + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[8][i] = 0x20212223 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[9][i] = 0x24252627 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[10][i] = 0x28292a2b + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[11][i] = 0x2c2d2e2f + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[12][i] = 0x30313233 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[13][i] = 0x34353637 + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[14][i] = 0x38393a3b + 0x10101010 * i;
		for (int i = 0; i < 8; i++)
			msg[15][i] = 0x3c3d3e3f + 0x10101010 * i;
		vsha256_transform(h, msg);
		for (int i = 0; i < 8; i++) {
			printf("lane %d: ", i);
			for (int j = 0; j < 8; j++)
				printf("%08x ", h[j][i]);
			printf("\n");
		}
	}
}
#endif
