// Simple dictionary implementation using open addresing, linear probing
// the input values are not hashed since we assume that keys have been
// already hashed (context: long message attack)

#include "dict.h"
#include "types.h"
//#include "util_char_arrays.h"
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



dict* dict_new(size_t nelements){
  /// initialize a dictionary that will hold nslots elements
  /// return a pointer to the dicitonary
  /// In memeory keys are seperate from the dictionary fields
  /// nslots in the dictionary will be alwais a power of two
  /// f is a pointer to a file that saves the 
  /// NOTE: key_size is fixed in as dict_key in types.h
  ///      It's 128 bit key

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
     d->slots[i].value = 0; // 0 if it is not occupied
     d->slots[i].key._uint64[0] = 0;
     d->slots[i].key._uint64[1] = 0;
  }  
   
   // printf("- Dictionary of size 0x%lx\n has been initialized\n", nslots);
  return d;
}


//void dict_add_element_to(dict* d, dict_key* key, size_t value, size_t key_size){
void dict_add_element_to(dict* d, uint64_t* key, size_t value){
  // add (key, value) to the dictionary, 0 <= value. Since we use value:=0 to indicate
  // that the associated slot in the dictionary is empty. In dictionary we store
  // value + 1
  // key is an array 64bit numbers
  // in this primitive design key is an array of bytes, thus key_size is how many bytes
  // value is a value that represents an index of some array thus it has the type size_t
  // dict = {slot1, slot2, ..., slot_n}

  /// our dictionary can be indexed using 64 bits
  // find which bin to put it in
  uint64_t h = key[0]; // for now we assume indices don't need more than 64 bits 
  h = h & (d->nslots - 1); // assumption nslots = 2^m; 
  // locate where to place (key, value)
  while (d->slots[h].value){ // value==0 means empty slot
    // linear probing
    // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    if (key[0] == d->slots[h].key._uint64[0]
	&& key[1] == d->slots[h].key._uint64[1]){
      // printf("a duplicated key has been detected at %lu\n", value);
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
  // Ensure the entered value is strictly greate than 0
  d->slots[h].value = value + 1;
  d->slots[h].key._uint64[0] = key[0];
  d->slots[h].key._uint64[1] = key[1];

}



size_t dict_get_value(dict* d, uint64_t* key){
  // we first need to find where to start probing

  size_t h =  key[0];
  h = h & (d->nslots - 1); // assumption nslots = 2^m <= 2^64; 

  // find (key, value) after 
  while (d->slots[h].value){// occupied slot
    // linear probing
    // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    if (key[0] == d->slots[h].key._uint64[0]
	&& key[1] == d->slots[h].key._uint64[1])
      return d->slots[h].value - 1;
    
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    h = (h+1) & (d->nslots - 1); // mod 2^nslots
  }

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  return 0; // no element is found
  
}


void dict_print(dict* d){
 

  for (size_t b=0; b<(d->nslots); ++b) {

    printf("slot=%lu, value=%lu, "
	   "key[1]key[0] = 0x%016lx%016lx\n",
	   b, d->slots[b].value,
	   d->slots[b].key._uint64[1], d->slots[b].key._uint64[0]);
    
    /* // print key */
    /* printf("key=0x"); */
    /* for (size_t k=0; k<key_size; ++k) */
    /*   printf("%x",(unsigned char) d->slots[b].key.bytes[k]); */
    /* printf("\n"); */
  }
  
}
 

int dict_has_key(dict* d, uint64_t* key){
  /// 1 if the dictionary has `key`
  /// 0 otherwise

  
 
  // we first need to find where to start probing
  size_t h =  key[0];
  h = h & (d->nslots - 1); // assumption nslots = 2^m; 

  // find (key, value) after 
  while (d->slots[h].value){ // value==0 means empty slot 
    // linear probing
     // current = element with index (h+i mod nslots)
    // puts("collision at adding element has been detected, linear probing");
    
    // check if key already exists
    
    if (key[0] == d->slots[h].key._uint64[0]
	&& key[1] == d->slots[h].key._uint64[1])
      // revert it to the original entered value
      return d->slots[h].value - 1;
    
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    h = (h+1) & (d->nslots - 1); // mod 2^nslots
  }

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  return 0; // no element is found
 
}

