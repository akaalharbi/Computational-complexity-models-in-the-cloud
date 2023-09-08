#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "../include/numbers_shorthands.h"
#include "../include/miller_rabin.h"
#include "../include/simple_dict.h"



#include <stdio.h>



// create dict
dict dict_create(u64 nelements)
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
void dict_insert(u64 val, dict* d)
{
  d->keys[val % d->p] = val;
}

// search
int dict_search(u64 val, dict* d)
{
  // if there is an element or not
  return (d->keys[val % d->p] > 0);
}


#ifdef TEST_MAIN
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
