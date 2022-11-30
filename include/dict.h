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
#include "util_char_arrays.h" // cmp_arrays
#include "shared.h"
#include "numbers_shorthands.h"

///-----------------------------------------------------///
///                  config                             ///
///-----------------------------------------------------///

// based on our experiments this will inlude > 91% of elements
// with filling rate 1. In conventional linear probing memory will
// be (1/0.91) * nslot * sizeof(uint32) while in our case it's nslots * sizeof(uint32)
#ifndef NPROBES_MAX
#define NPROBES_MAX 32
#endif

#ifndef FILLING_RATE
#define FILLING_RATE 1
#endif



///-----------------------------------------------------///
///                  data structure                     ///
///-----------------------------------------------------///
typedef struct  __attribute__((aligned(32))) {
  uint64_t* keys   __attribute__((aligned(32)));
  size_t nslots; // number of elements  in the dictionary
  size_t nbuckets; // a bucket contains x0 slots depending avx register length
  size_t nslots_per_bucket; // = nslots / nbuckets
  //size_t scale; // = log2(nslots_per_bucket)
  size_t nprobes_insert;
  size_t nprobes_lookup;
  size_t nelments_succ_lookup;
} dict;




//-----------------------------------------------------//
//                    functions                        //
//-----------------------------------------------------//

dict *dict_new(size_t nslots);
void dict_init(dict* d);
void dict_free(dict* d);
size_t  dict_memory(size_t nelements);

int dict_add_element_to(dict* d, uint64_t idx, u32 val);
u32 dict_get_value(dict *d, uint64_t store_as_idx, u32 val);
void dict_print(dict *d);




#endif



