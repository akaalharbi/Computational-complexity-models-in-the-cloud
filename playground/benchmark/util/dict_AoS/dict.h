#ifndef DICT_LINEAR_PROBING
#define DICT_LINEAR_PROBING


#include "types.h"
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


typedef struct {
  /// Although our ideal case is to store 96-bits digest as a keys
  /// However, we prefer to save the space and consider only 60 bits
  /// the long message attack will deal with the false positives
  
  uint64_t key; // the key has a fixed length given at the beginning of the run
  size_t value; // if value==0 then its an empty slot. Any element added to dict
                // need to be incremented by one.
  // no need for this
  // int is_occupied; // 1 bit. 0 if it is empty, 1 if slot has an element
  
} slot;
 

typedef struct {
  size_t nslots; // number of elements  in the dictionary
  size_t nprobes_insert;
  size_t nprobes_lookup;
  slot* slots; // == slot slots[nslots] but on the heap
  
} dict;






//-----------------------------------------------------//
//                    functions                        //
//-----------------------------------------------------//
slot *slot_new();
dict *dict_new(size_t nslots);
void dict_free(dict *d);
int dict_memory(size_t nelements);

void dict_add_element_to(dict* d, uint64_t key[2], size_t value);
size_t dict_get_value(dict* d, uint64_t key[2]);
void dict_print(dict *d);

#endif
