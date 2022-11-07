// SoA dictionary
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
#define ALIGNMENT 32
#include <immintrin.h>


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

  //dict* d = (dict *) malloc(sizeof(dict));
  dict* d = (dict *) aligned_alloc(ALIGNMENT, (sizeof(dict)));
  // load factor = 1/2 =>
  // E[#probing] = 1.5 hits
  // E[#probing] = 2.5 misses
  
  // Use a defined filling rate type.h, empiricall data shows 0.77
  // is the best
  size_t nslots = (size_t)  ceil((1/FILLING_RATE)*nelements);
  
  // the number of slots is multiple of alignment
  // this is for SIMD. Don't confuse it with sha256 padding
  // multiple of ALIGNMENT BYTES
  int padding_alignment =   nslots %  (ALIGNMENT/sizeof(uint64_t)) ;
  padding_alignment = (ALIGNMENT/sizeof(uint64_t)) - padding_alignment ;
  //nslots = nslots + padding_alignment;
  d->nslots = nslots;
  printf("nslots=%lu, padding=%d, alignmen/sizeof(uint64)=%lu\n",
	 d->nslots, padding_alignment, (ALIGNMENT/sizeof(uint64_t)));

  d->nprobes_insert=0;
  d->nprobes_lookup=0;
  d->keys = (uint64_t*) aligned_alloc(ALIGNMENT,
		         (padding_alignment + nslots)*(sizeof(uint64_t)));
  // (uint64_t*) malloc(d->nslots*(sizeof(uint64_t)));
  d->values = (uint64_t*) aligned_alloc(ALIGNMENT,
			  (padding_alignment + nslots)*(sizeof(uint64_t)));
				
  d->memory_estimates = d->nslots*(sizeof(uint64_t)) // keys
                      + d->nslots*(sizeof(uint64_t)) // values
                      + sizeof(dict);
  
  // reserve space in memeory
  // Ensure the keys are next to each other
   for (size_t i = 0; i < nslots; ++i) {
     d->keys[i] = 0; // 0 if it is not occupied
     d->values[i] = 0;
  }  
   // printf("- Dictionary of size 0x%lx\n has been initialized\n", nslots);
  return d;
}

void dict_free(dict* d){
  free(d->keys);
  free(d->values);
  //free(d);
}


int dict_memory(size_t nelements){
  /// return memory estimation of the dictionary size
  int estimate = (size_t)  ceil((1/FILLING_RATE)*nelements);
  int padding_alignment = (ALIGNMENT/sizeof(uint64_t)) - nelements%( (ALIGNMENT/sizeof(uint64_t)) );
  // how many slots in the dictionary
  estimate = estimate + padding_alignment;

  estimate = estimate*(sizeof(uint64_t)) // keys
                      + estimate*(sizeof(uint64_t)) // values
                      + sizeof(dict);
   return estimate;
}
//void dict_add_element_to(dict* d, dict_key* key, size_t value, size_t key_size){
void dict_add_element_to(dict* d, uint64_t key[2], size_t value){
  // add (key, value) to the dictionary, 0 <= value. Since we use value:=0 to indicate
  // that the associated slot in the dictionary is empty. In dictionary we store
  // value + 1
  // key is an array 64bit numbers
  // in this primitive design key is an array of bytes, thus key_size is how many bytes
  // value is a value that represents an index of some array thus it has the type size_t
  // dict = {slot1, slot2, ..., slot_n}

  /// our dictionary can be indexed using 64 bits
  // find which bin to put it in
  // apologies: only store 64 bit 
  uint64_t h = key[0]; // for now we assume indices don't need more than 64 bits 
  h = h % d->nslots; // assumption nslots = 2^m; 
  // locate where to place (key, value)
  while (d->values[h]){ // value==0 means empty slot
    // linear probing
    //  ++i;
    // current = &d->slots[  (h+i) & (d->nslots - 1)  ];
    // reference: The art of simd programming, Sergey Slotin
    ++h;
    if (h >= d->nslots)
      h = h - d->nslots;
    
    
    #ifdef NPROBES_COUNT
    ++(d->nprobes_insert);
    #endif
  }


  /// Found an empty slot, update its 
  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  // Ensure the entered value is strictly greate than 0
  d->values[h] = value + 1;
  d->keys[h] = key[0];
}



void print_m25i(__m256i a, char* text){
  uint64_t A[4] = {0};
  _mm256_storeu_si256((__m256i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<4; ++i) {
    printf("%016lx, ", A[i]);
  }
  puts("");
}



size_t dict_get_value(dict* d, uint64_t key[2]){
  // we first need to find where to start probing
  // apologies: only check 64
  //int step = ALIGNMENT/sizeof(d->keys[0]); // for the while loop
  int step = ALIGNMENT/sizeof(uint64_t); // for the while loop
  size_t h =  key[0];
  // h = h & (d->nslots - 1); // assumption nslots = 2^m <= 2^64; 
  h = h % d->nslots;
  // to access an element with memory address that is multiple of ALIGNMENT
  //h = h - (h%step);
  h = h - (h&(step-1)); // since step = 2^r for some r 
  // min(h, nslots - alignement) so we can get #alignement*bytes 
  // h = (h > d->nslots - ALIGNMENT) ? d->nslots - ALIGNMENT : h;
  // no need for this 

  int is_key_found = 0;
  int has_empty_slot = 0; //1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);
  
  __m256i dict_keys_simd;// = _mm256_loadu_si256((__m256i*)  &(d->keys[h]));
  __m256i lookup_key_simd = _mm256_setr_epi64x(key[0], key[0], key[0], key[0]);
  __m256i zero_vect = _mm256_setzero_si256();
  __m256i comp_vect_simd;

  
    
  
  // find (key, value) after
  // 1, 1, 0, 1
  // 1, 1, 0, 1
  
  while (!has_empty_slot){// occupied slot,
    // we are relying that the key != 0, this is a negligible event
    // linear probing

    // get new fresh keys
    //printf("keys=%p\n", &d->keys[h]);
    
    dict_keys_simd = _mm256_load_si256((__m256i*)  &(d->keys[h]));
      /// -----------------------------------------------///
     ///                   TEST 1                       ///
    ///------------------------------------------------///
    /* does key equal one of the slots */
    comp_vect_simd = _mm256_cmpeq_epi64(lookup_key_simd, dict_keys_simd);
    is_key_found = 1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);

    #ifdef VERBOSE_LEVEL
    printf("step=%d, h=%lu, nslots=%lu\n", step, h, d->nslots);
    print_m25i(dict_keys_simd, "dict_keys_simd");
    print_m25i(lookup_key_simd, "lookup_key_simd");
    print_m25i(comp_vect_simd, "compare with keys");
    printf("is_key_found? %d\n", is_key_found);
    #endif
    
    // todo start here                                  
    if (is_key_found) {
      // movemask_ps(a) will return 4 bytes m, where m[i] := some mask with ith 64 bit entry
      int loc = _mm256_movemask_epi8(comp_vect_simd);
      // some magical formula, too tired to explain it's 00:12am
      size_t idx_val = (__builtin_ctz(loc)>>3) + h;
      // since we deal with values of size 64 will be repeated twice
      return d->values[idx_val];// - 1; // let caller deal with correcting the offset
    }

    ///                   TEST 2                           ///
    ///----------------------------------------------------///
    /// is one the slots empty? then no point of probing
    comp_vect_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);
    // _mm256_testz_si256 will return if comp_vect_simd AND comp_vect_simd = 0 as int
    has_empty_slot = 1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);




    // update the index for keys load
    h += step;
    if (h >= d->nslots)
      h = 0;

    #ifdef VERBOSE_LEVEL
    print_m25i(dict_keys_simd, "dict_keys_simd");
    print_m25i(zero_vect, "zero_vect");
    print_m25i(comp_vect_simd, "compare with zeros");
    printf("has empty slot? %d\n", has_empty_slot);
    #endif


    #ifdef VERBOSE_LEVEL
    print_m25i(dict_keys_simd, "fresh dict_keys_simd");
    printf("has empty slot? %d\n", has_empty_slot);
    #endif
    
    #ifdef NPROBES_COUNT
    ++(d->nprobes_lookup);
    #endif
    //printf("inside next: step=%d, h=%lu, nslots=%lu\n", step, h, d->nslots);
  }

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
  return 0; // no element is found
  
}


void dict_print(dict* d){
 

  for (size_t b=0; b<(d->nslots); ++b) {

    printf("slot=%lu, value=%lu, "
	   "key = 0x%016lx\n",
	   b, d->values[b],
	   d->keys[b]);
    
    /* // print key */
    /* printf("key=0x"); */
    /* for (size_t k=0; k<key_size; ++k) */
    /*   printf("%x",(unsigned char) d->slots[b].key.bytes[k]); */
    /* printf("\n"); */
  }
  
}
 
/// legacy code 
/* int dict_has_key(dict* d, uint64_t key[2]){ */
/*   /// 1 if the dictionary has `key` */
/*   /// 0 otherwise */

  
 
/*   // we first need to find where to start probing */
/*   size_t h =  key[0]; */
/*   // h = h & (d->nslots - 1); // assumption nslots = 2^m;  */
/*   h = h % d->nslots; */
/*   // find (key, value) after  */
/*   while (d->slots[h].value){ // value==0 means empty slot  */
/*     // linear probing */
/*      // current = element with index (h+i mod nslots) */
/*     // puts("collision at adding element has been detected, linear probing"); */
    
/*     // check if key already exists */
    
/*     if (key[0] == d->slots[h].key */
/* 	&& key[1] == d->slots[h].key._uint64[1]) */
/*       // revert it to the original entered value */
/*       return d->slots[h].value - 1; */
    
/*     // current = &d->slots[  (h+i) & (d->nslots - 1)  ]; */
/*     // h = (h+1) & (d->nslots - 1); // mod 2^nslots */
/*         // current = &d->slots[  (h+i) & (d->nslots - 1)  ]; */
/*     ++h; */
/*     if (h >= d->nslots) */
/*       h = h - d->nslots; */

/*   } */

/*   // update current->key = key  */
/*   // memcpy(current->key, key->bytes, key_size); */
/*   return 0; // no element is found */
 
/* } */

