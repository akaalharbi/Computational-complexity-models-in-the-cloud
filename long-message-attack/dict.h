#ifndef DICT_LINEAR_PROBING
#define DICT_LINEAR_PROBING


#include "types.h"
// #include "sha256.h" // digest 
#include <stddef.h>  // size_t
#include <stdint.h> // uintXY_t
#include <stdlib.h>  // malloc
#include <string.h> // memcpy
#include <stdio.h>  // printf
#include "types.h" // digest union type 
#include "util_char_arrays.h" // cmp_arrays

///-----------------------------------------------------///
///                  data structure                     ///
///-----------------------------------------------------///

typedef struct {
  /// a slot can hold one element at most!
  /// This implementaion suggest that slot1's key and slot2's key are
  /// NOT next to each other in the memory. Durint the dictionary init
  /// we will ensure that they are contagious in the memory. 

  
  int is_occupied: 1; // 1 bit. 0 if it is empty, 1 if slot has an element
  digest* key; // the key has a fixed length given at the beginning of the run
  size_t value; // to be edited when we use a digest that has more than 
} slot;


typedef struct {
  UINT nslots; // number of elements  in the dictionary
  slot* slots; // == slot slots[nslots] but on the heap
} dict;






//-----------------------------------------------------//
//                    functions                        //
//-----------------------------------------------------//
slot *slot_new();
dict *dict_new(size_t nslots, size_t key_size);
void dict_add_element_to(dict* d, digest* key,size_t value, size_t key_size);
int dict_has_key(dict* d, digest* key, size_t key_size);
size_t dict_get_value(dict* d, digest* key, size_t key_size);
void dict_print(dict *d, size_t key_size);

#endif
