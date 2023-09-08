
#ifndef NPROBES_MAX
#define NPROBES_MAX 32
#endif

#ifndef FILLING_RATE
#define FILLING_RATE 1
#endif



#ifndef SIMPLE_DICT
#define SIMPLE_DICT



#include "numbers_shorthands.h"
#include <stddef.h>

typedef struct {
  /* simple dictionary only has an array and a prime that p=|dict| */
  u64* keys;
  u64 nslots;
  u64 nelements;
  u64 nelements_asked_to_be_inserted;
} dict;


// create dict
dict dict_new(u64 nelements);
// insert
void dict_add_element_to(dict* d, u8* state);
// search
int dict_has_elm(dict* d, u8* state);
void dict_free(dict *d);
size_t  dict_memory(size_t nelements);


#endif
