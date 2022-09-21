// Simple dictionary implementation using open addresing, linear probing
// the input values are not hashed since we assume that keys have been
// already hashed (context: long message attack)

#include "dict.h"
#include "sha256.h"
#include "types.h"
#include "util_char_arrays.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"



/* //-----------------------------------------------------// */
/* //                  data structure                     // */
/* //-----------------------------------------------------// */
/// See the header file

//-----------------------------------------------------//
//                       methods                       //
//-----------------------------------------------------//



dict* dict_new(size_t nelements, size_t key_size){
  /// initialize a dictionary that will hold nslots elements
  /// return a pointer to the dicitonary
  /// In memeory keys are seperate from the dictionary fields
  /// nslots in the dictionary will be alwais a power of two

  dict* d = (dict *) malloc(sizeof(dict));
  // load factor = 1/2 =>
  // E[#probing] = 1.5 hits
  // E[#probing] = 2.5 misses
  
  // also guarantees that nslots is a power of 2 
  size_t nslots = (size_t) 1 << ( (size_t) ceil(log2( nelements) + 1 ) );
  d->nslots = nslots;
  
  // reserve space in memeory
  d->slots = malloc(sizeof(slot)*nslots);
  // Ensure the keys are next to each other
   for (size_t i = 0; i < nslots; ++i) {

    // each key reserves `key_size` bytes, thus we need to move key_size steps
    // to go to the next slot key
     //d->slots[i].key.bytes;
     d->slots[i].is_occupied = 0; // to be removed 
     d->slots[i].value = 0; // 0 if it is not occupied
  }

  printf("- Dictionary of size 0x%lx\n has been initialized\n", nslots);
  return d;
}


void dict_add_element_to(dict* d, dict_key* key, size_t value, size_t key_size){
  // todo edit this 
  // add (key, value) to the dictionary
  // in this primitive design key is an array of bytes, thus key_size is how many bytes
  // value is a value that represents an index of some array thus it has the type size_t
  // todo if we entered a key twice, there will be two entries :(

  // dict = {slot1, slot2, ..., slot_n}
  // find which bin to put it in
  uint64_t h = key->_uint64[0]; // for now we assume indices don't need more than 64 bits 
  h = h & (d->nslots - 1); // assumption nslots = 2^m; 
  // locate where to place (key, value)
  while (d->slots[h].is_occupied){
    // linear probing
    // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    if (key->_uint64[0] == d->slots[h].key._uint64[0]){
      printf("a duplicated key has been detected at %lu\n", value);
      is_there_duplicate = 1;
      return;
    }
    //  ++i;
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    h = (h+1) & (d->nslots - 1);
  }

  /// Found an empty slot, update its 
  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  d->slots[h].is_occupied = 1;
  d->slots[h].value = value;
  d->slots[h].key._uint64[0] = key->_uint64[0];

}



size_t dict_get_value(dict* d, dict_key* key, size_t key_size){
  // we first need to find where to start probing

  size_t h =  key->_uint64[0];
  h = h & (d->nslots - 1); // assumption nslots = 2^m; 

  // find (key, value) after 
  while (d->slots[h].is_occupied){
    // linear probing
    // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    if (key->_uint64[0] == d->slots[h].key._uint64[0])
      return d->slots[h].value;
    
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    h = (h+1) & (d->nslots - 1); // mod 2^nslots
  }

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  return 0; // no element is found
  
}


void dict_print(dict* d, size_t key_size){


  for (size_t b=0; b<(d->nslots); ++b) {

    printf("slot=%lu, value=%lu, occupied=%d, ", b, d->slots[b].value, d->slots[b].is_occupied);
    // print key
    printf("key=0x");
    for (size_t k=0; k<key_size; ++k)
      printf("%x",(unsigned char) d->slots[b].key.bytes[k]);
    printf("\n");
  }
  
}
 

int dict_has_key(dict* d, dict_key* key, size_t key_size){
  /// 1 if the dictionary has `key`
  /// 0 otherwise

  

  // we first need to find where to start probing
  size_t h =  key->_uint64[0];
  h = h & (d->nslots - 1); // assumption nslots = 2^m; 

  // find (key, value) after 
  while (d->slots[h].is_occupied){
    // linear probing
    // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    if (key->_uint64[0] == d->slots[h].key._uint64[0])
      return d->slots[h].value;
    
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    h = (h+1) & (d->nslots - 1); // mod 2^nslots
  }

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  return 0; // no element is found
 
}

