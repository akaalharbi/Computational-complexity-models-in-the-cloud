// SoA dictionary
// Simple dictionary implementation using open addresing, linear probing
// the input values are not hashed since we assume that keys have been
// already hashed (context: long message attack)

#include "dict.h"
#include "config.h"
//#include "util_char_arrays.h"
#include <emmintrin.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"

#include <immintrin.h>

//  how many bits we store as a value in dictionary
#define NBITS_VAL 32 // To avoid having magical numbers in the code



/* //-----------------------------------------------------// */
/* //                  data structure                     // */
/* //-----------------------------------------------------// */
/// See the header file

//-----------------------------------------------------//
//                       methods                       //
//-----------------------------------------------------//
/* represent n in <= 6 char  */

dict* dict_new(size_t nelements){
  /// dict v4: we take part of the value as in index and store the rest
  ///          as a value. The value is maximally 32 bits 
  /// dict v3: we don't store values, we only store 64bit of the key
  ///         in d->keys[k]. When d->keys[k]=0, it means that we have
  ///         we have an empty slot


  dict* d = (dict *) aligned_alloc(ALIGNMENT, (sizeof(dict)));

  // on my laptop 8 elements per bucket
  int nslots_per_bucket = AVX_SIZE / NBITS_VAL; // we store 32 bits per value
  //size_t scale = (int) log2(nslots_per_bucket);
  // we will add few slots to enusre 
  // size_t nslots = (size_t)  ceil((1/FILLING_RATE)*nelements);
  size_t nslots = nelements;

  // nsolts = nbuckets * nslots_per_bucket; nbuckets divides nslots
  nslots = nslots + (-nslots % nslots_per_bucket);
  //+ todo does negative gives the correct result here?

  /// Save configured variables as dictionary entries 
  d->nslots = nslots;
  d->nbuckets = nslots/nslots_per_bucket;
  d->nslots_per_bucket = nslots_per_bucket;
  //d->scale = scale;
  d->nprobes_insert=0;
  d->nprobes_lookup=0;
  d->nelments_succ_lookup = 0;
  d->keys = (uint64_t*) aligned_alloc(ALIGNMENT,
				      (nslots)*(sizeof(uint64_t)));



  // Ensure keys are zeroed
  #pragma omp simd
  for (size_t i = 0; i < nslots; ++i) 
     d->keys[i] = 0; // 0 if it is not occupied

  return d;
}

inline void dict_free(dict* d){
  free(d->keys);
}


size_t dict_memory(size_t nelements){
  /// return memory estimation of the dictionary size
  int nslots_per_bucket = AVX_SIZE / NBITS_VAL; // we store 32 bits per value
  size_t nslots = nelements;
  nslots = nslots + (-nslots % nslots_per_bucket);
  size_t estimate = nslots*(sizeof(uint32_t)) + sizeof(dict);
  
  return estimate/1000.0;
}


int dict_add_element_to(dict* d, uint64_t store_as_idx, uint32_t val){
  // =========================================================================+
  // returns 1 if an element has been added, 0 otherwise                      |
  // This dictionary is unusual:                                              |
  // User have Value = (64bit) || (32bit), the user choses to split as they   |
  // wish. The reason, we don't do split by ourselves is that we might have   |
  // clustering if gcd(nbuckets, store_as_idx) > 1, this is likely to happen  |
  // due to the fact that store_as_idx = 2^b * g.                             |
  // In the future we might write a cleaner interface that only takes value   |
  // and dict*                                                                |
  // -------------------------------------------------------------------------|
  // Note about buckets: 
  //==========================================================================+

  /// Use linear probing to add element to the array d->keys
  // get bucket number, recall keys[nbuckets*nslots_per_bucket]
  uint64_t idx = store_as_idx % d->nbuckets; 
  // todo aestheitcs: which version do you like more?
  // idx = idx << (d->scale) ; // get the index of bucket in keys array
  idx = idx * (d->nslots_per_bucket) ; // get the index of bucket in keys array
  // linear probing 
  // locate where to place (key, value)
  // in insertion we are not going to use buckets explicitly
  for (int i=0; i<NPROBES_MAX; ++i) {
    // found an empty slot inside a bucket
    if (d->keys[idx] == 0) { // found an empty slot
      d->keys[idx] = val;
      return 1;
    }


    ++idx;
    // reduce mod n->slots
    if (idx>d->nslots)
      idx = 0;
  }
  return 0; // element has been added
}


/* This method should be sent to some other file */
void print_m25i(__m256i a, char* text){
  uint64_t A[4] = {0};
  _mm256_storeu_si256((__m256i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<4; ++i) {
    printf("%016lx, ", A[i]);
  }
  puts("");
}



uint32_t dict_get_value(dict *d, uint64_t store_as_idx, uint32_t val){
  // =========================================================================+
  // This dictionary is unusual:                                              |
  // User have Value = (64bit) || (32bit), the user choses to split as they   |
  // wish. The reason, we don't do split by ourselves is that we might have   |
  // clustering if gcd(nbuckets, store_as_idx) > 1, this is likely to happen  |
  // due to the fact that store_as_idx = 2^b * g.                             |
  // In the future we might write a cleaner interface that only takes value   |
  // and dict*                                                                |
  //==========================================================================+
  
  uint64_t h = store_as_idx; // for now we assume indices don't need more than 64 bits 
  h = h % d->nbuckets; 
  size_t idx = h * d->nslots_per_bucket ; // get the index of bucket in keys array
  


  int is_key_found = 0;
  // it's enough to check if the first element is empty
  // in this version we don't need to check if we hit zero or not.
  /* int empty_bucket = 0; //1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd); */
  // we can remove one of the above variables 
  
  __m256i dict_keys_simd;// = _mm256_loadu_si256((__m256i*)  &(d->keys[h]));
  __m256i lookup_key_simd = _mm256_set1_epi32(val); // (val, val, ..., val) 8times
  //__m256i zero_vect = _mm256_setzero_si256(); // no need for this with buckets
  __m256i comp_vect_simd;

  
  
  
  // loop at most NPROBES_MAX/SIMD_LEN since we load SIMD_LEN
  // elements from dictionary each loop.
  for (size_t i=0; i< (int) (NPROBES_MAX/SIMD_LEN); ++i) {
        
    // we are relying that the val != 0, Pr(val !=0 ) = 1/(2^32)
    // linear probing

    // get new fresh keys from one bucket
    dict_keys_simd = _mm256_load_si256((__m256i*)  &(d->keys[idx]));
    // -----------------------------------------------//
    //                   TEST 1                       //
    /*  Does key equal one of the slots?              */
    //------------------------------------------------//
    comp_vect_simd = _mm256_cmpeq_epi32(lookup_key_simd, dict_keys_simd);
    is_key_found = 1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);

    #ifdef VERBOSE_LEVEL
    printf("step=%d, h=%lu, nslots=%lu\n", step, h, d->nslots);
    print_m25i(dict_keys_simd, "dict_keys_simd");
    print_m25i(lookup_key_simd, "lookup_key_simd");
    print_m25i(comp_vect_simd, "compare with keys");
    printf("is_key_found? %d\n", is_key_found);
    #endif
    

    if (is_key_found) { // @todo explains theese magical lines!
      // movemask_ps(a) will return 4 bytes m, where m[i] := some mask with ith 64 bit entry

      // wi:32-bit,  w8 w7 w6 w5 w4 w3 w2 w1 - move mask ->
      // (hw8||hw7) || (hw6||hw5) || (hw4||hw3) || (hw2||hw1) 
      // hwi: 4-bits extracted as MSB of each consecutive 8bits in wi
      uint32_t loc = _mm256_movemask_epi8(comp_vect_simd);

      // We wish to find the wi that is not 0.
      // each wi on the right that equals zero, will add 4 zeros
      // (since movemask gets 4 bits from each 32-bit word)
      // __builtin_ctz(x): Returns the number of trailing 0-bits in x, starting
      // at the least significant bit position. undefined if x==0.
      size_t idx_val = (__builtin_ctz(loc)>>2) + idx;
      // note in our case loc != 0 since a key was found.
      return d->keys[idx_val];
    }

    // -------------------------------------------------//
    //                   TEST 2                         //
    /* is one the slots empty? then no point of probing */
    //- ------------------------------------------------//
    // no need for the commented instruction since we only check the first value
    // comp_vect_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);
    // copy the value of the first slot of the bucket, check is it 0?

    // with this new version we don't need to check if the first element is 0
    // for two reasons: 1- The loop size is already defined.
    // 2- if we hit zero but kept moving then worst case scenario we wil have a false positive
    // we just need to increase the number of needed hits.
    
    /* empty_bucket = ( _mm256_cvtsi256_si32(comp_vect_simd) == 0 ); */

    /* if (empty_bucket) */
    /*   return 0; */


    // Linear probing
    // update the index for keys load
    idx += d->nslots_per_bucket; // move to the next bucket
    if (idx >= d->nslots)
      idx = 0;

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


  return 0; // no element is found
  
}



void dict_print(dict* d){
 
  for (size_t b=0; b<(d->nslots); ++b) {
    printf("slot=%lu, "
	   "key = 0x%016lx\n",
	   b,
	   d->keys[b]);
  }
}
 
