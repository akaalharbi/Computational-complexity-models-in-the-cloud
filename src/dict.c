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
  size_t nslots = (size_t)  ceil((1/FILLING_RATE)*nelements);
  

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
  size_t nslots = (size_t)  ceil((1/FILLING_RATE)*nelements);
  nslots = nslots + (-nslots % nslots_per_bucket);
  size_t estimate = nslots*(sizeof(uint32_t)) + sizeof(dict);
  
  return estimate/1000.0;
}


void dict_add_element_to(dict* d, uint64_t store_as_idx, uint32_t val){
  // =========================================================================+
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
      return;
    }


    ++idx;
    // reduce mod n->slots
    if (idx>d->nslots)
      idx = 0;
  }
  // done
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
  int empty_bucket = 0; //1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);
  __m256i dict_keys_simd;// = _mm256_loadu_si256((__m256i*)  &(d->keys[h]));
  __m256i lookup_key_simd = _mm256_set1_epi32(val); // (val, val, ..., val) 8times
  __m256i zero_vect = _mm256_setzero_si256();
  __m256i comp_vect_simd;

  
  
  // while (!empty_bucket){// has occupied slot,
  for (size_t i=0; i<NPROBES_MAX; ++i) {
        
    // we are relying that the val != 0, Pr(val !=0 ) = 1/(2^32)
    // linear probing

    // get new fresh keys from one bucket
    dict_keys_simd = _mm256_load_si256((__m256i*)  &(d->keys[idx]));
    // -----------------------------------------------//
    //                   TEST 1                       //
    /*  Does key equal one of the slots?              */
    //------------------------------------------------//
    comp_vect_simd = _mm256_cmpeq_epi64(lookup_key_simd, dict_keys_simd);
    is_key_found = 1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd);

    #ifdef VERBOSE_LEVEL
    printf("step=%d, h=%lu, nslots=%lu\n", step, h, d->nslots);
    print_m25i(dict_keys_simd, "dict_keys_simd");
    print_m25i(lookup_key_simd, "lookup_key_simd");
    print_m25i(comp_vect_simd, "compare with keys");
    printf("is_key_found? %d\n", is_key_found);
    #endif
    
    //+ todo: correct this area                                  
    if (is_key_found) { // @todo explains theese magical lines!
      // movemask_ps(a) will return 4 bytes m, where m[i] := some mask with ith 64 bit entry
      int loc = _mm256_movemask_epi8(comp_vect_simd);
      // some magical formula, too tired to explain it's 00:12am
      size_t idx_val = (__builtin_ctz(loc)>>3) + h;
      // since we deal with values of size 64 will be repeated twice
      return d->keys[idx_val];// - 1; // let caller deal with correcting the offset
    }

    // -------------------------------------------------//
    //                   TEST 2                         //
    /* is one the slots empty? then no point of probing */
    //- ------------------------------------------------//
    
    comp_vect_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);
    empty_bucket = _mm256_movemask_epi8(comp_vect_simd);

    if (empty_bucket)
      return 0;



    // update the index for keys load
    //!
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

  // update current->key = key 
  // memcpy(current->key, key->bytes, key_size);
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
 
