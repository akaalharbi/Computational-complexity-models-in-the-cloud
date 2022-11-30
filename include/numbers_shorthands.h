#ifndef NEAT_NUMBER_SHORTHANDS
#define NEAT_NUMBER_SHORTHANDS
#include <stdint.h>

typedef unsigned char          u8  ;
typedef unsigned short	       u16 ;
typedef unsigned int	       u32 ;
typedef unsigned long long     u64 ;

#ifdef __SIZEOF_INT128__
typedef unsigned __int128      u128;
#endif

typedef signed char	        s8 ;
typedef short			s16;
typedef int			s32;
typedef long long		s64;



#endif
