// SoA, the keys will have an indenependent array from 
// values array
#ifndef DICT_LINEAR_PROBING
#define DICT_LINEAR_PROBING


#include "config.h"
// #include "sha256.h" // dict_key 
#include <stddef.h>  // size_t
#include <stdint.h> // uintXY_t
#include <stdlib.h>  // malloc
#include <string.h> // memcpy
#include <stdio.h>  // printf
#include "types.h" // dict_key union type 
#include "util/util_char_arrays.h" // cmp_arrays
#include "shared.h"
///-----------------------------------------------------///
///                  data structure                     ///
///-----------------------------------------------------///


/// legacy code 
/* typedef union { */
/*   /// the keys have fixed size */
/*   /// Wrap state type in a union for easier handling outside the sha256 function  */
/*   char bytes[16]; */
/*   unsigned int _uint32[4]; */
/*   uint64_t _uint64[2];  */

/* }  dict_key; // 128 bit */



 

typedef struct  __attribute__((aligned(32))) {
  uint64_t* keys   __attribute__((aligned(32)));
  size_t nslots; // number of elements  in the dictionary
  size_t nprobes_insert;
  size_t nprobes_lookup;
  //  size_t memory_estimates;
} dict;




//-----------------------------------------------------//
//                    functions                        //
//-----------------------------------------------------//

dict *dict_new(size_t nslots);
void dict_init(dict* d);
void dict_free(dict* d);
size_t  dict_memory(size_t nelements);
void dict_add_element_to(dict* d, uint64_t key[2]);
#pragma omp declare simd uniform( d )
void dict_get_values_simd(dict* d, uint64_t keys[4], uint64_t found_keys[4]);
void dict_print(dict *d);

#ifndef FILLING_RATE
#define FILLING_RATE 0.77
#endif

#endif



