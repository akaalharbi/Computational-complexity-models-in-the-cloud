#ifndef SIMPLE_DICT
#define SIMPLE_DICT

#include "numbers_shorthands.h"

typedef struct {
  /* simple dictionary only has an array and a prime that p=|dict| */
  u64* keys;
  u64 p;
} dict;

dict dict_create(u64 nelements);
void dict_insert(u64 val, dict* d);
int dict_search(u64 val, dict* d);


#endif
