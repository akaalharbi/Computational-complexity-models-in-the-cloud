#ifndef VSHA256
#define VSHA256


typedef unsigned char		u8;
typedef unsigned short		u16;
typedef unsigned int		u32;
typedef unsigned long long	u64;
typedef signed char		s8;
typedef short			s16;
typedef int			s32;
typedef long long		s64;

/*****************************************************************************/

/*
 * This is the AVX-2 specific part. It defines a "vector" type
 */

typedef u32 vector __attribute__((vector_size(32)));	// 32-byte vectors (AVX2)
#define VECTOR_ALIGNED  __attribute__ ((aligned (32)))	// AVX-2 vectors must be aligned
#define AVX_ALIGNED VECTOR_ALIGNED
/*****************************************************************************/


/* functions */
vector vbroadcast(u32 x);
void __sha256_init(vector state[8]);
void vsha256_setup();
void __sha256_transform(vector state[], const vector W[]);
void vsha256_init(u32 state[8][8]);
void vsha256_transform(u32 state[8][8], const u32 W[16][8]);
#endif
