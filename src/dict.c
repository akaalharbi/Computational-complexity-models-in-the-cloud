#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "../include/numbers_shorthands.h"
#include "../include/miller_rabin.h"
#include "../include/dict.h"



#include <stdio.h>



// create dict
dict dict_new(u64 nelements)
{
  u64 p = 1;
  /* Find the largest prime smaller than nelements */
  for (u64 i=nelements; i>0; --i) {
    if (Miller(i, 10)){
      p = i;
      i = 0;
      break;
    }
  }

  // allocate array:
  size_t nbytes = sizeof(u64)*p;
  u64* keys = (u64*) malloc(nbytes);
  // set all elements to 0;
  memset(keys, 0, nbytes);
  dict d = {keys, p};
  printf("Initialized a dict with nelements=%llu\n", p);
  
  return d;
}



// insert
void dict_add_element_to(dict* d, u8* state)
{
  /* take the first 8 bytes */
  u64 val = ((u64*) state)[0];
  // only increase the number of elements if the slot is empty
  if (d->keys[val % d->nslots] == 0) {d->nelements += 1;}
  // write the value on the slot even if it has a previous value
  d->keys[val % d->nslots] = val;
  d->nelements_asked_to_be_inserted += 1;
}

// search
int dict_has_elm(dict* d, u8* state)
{

  // a code is not complete without this gymnastic
  u64 val = ((u64*) state)[0];
  // if there is an element or not
  return (d->keys[val % d->nslots] > 0);
}



void dict_free(dict *d) { free(d->keys); }

size_t  dict_memory(size_t nelements){ return sizeof(u64)*nelements ;}


void dict_print(dict *d){return;}


#ifdef TEST_MAIN_SIMPLE_DICT
#include "../include/test_dict.h"

int main()
{
  u64 n = 40000;
  dict d = dict_create(32*n);
  u64 nfound = 0;
  
  for (u64 i = 0; i<n; ++i)
    dict_insert(L_insert[i], &d);

  for (u64 i = 0; i<n; ++i){
    if (dict_search(L_search[i], &d))
      nfound += 1;
    
  }

  printf("nfound=%llu, i.e. %0.2f%%\n", nfound,100*nfound/(20055.0));

}
#endif
