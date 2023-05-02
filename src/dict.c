// SoA dictionary
// Simple dictionary implementation using open addresing, linear probing
// the input values are not hashed since we assume that keys have been
// already hashed (context: long message attack)

#include "confg_math_func.h"
#include "numbers_shorthands.h"
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
#include "util_char_arrays.h"
#include <sys/mman.h> 
#include <immintrin.h>

//  how many bits we store as a value in dictionary


void print_m256i(__m256i a, char* text){
  uint32_t A[8] = {0};
  _mm256_storeu_si256((__m256i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<8; ++i) {
    printf("%02x, ", A[i]);
  }
  puts("");
}



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

  int nslots_per_bucket = SIMD_LEN; /* See ../include/config.h */
  size_t nslots = nelements;

  //+ todo does negative gives the correct result here?
  //nslots = nslots + (-nslots % nslots_per_bucket);


  /// Save configured variables as dictionary entries 

  d->nbuckets = (nelements/nslots_per_bucket) + 1;
  d->nslots_per_bucket = nslots_per_bucket;
  d->nslots = (d->nbuckets)*(d->nslots_per_bucket);

  /* for huge pages : required memory for nslots is a mutlipe GPAGE_SIZE  */
  /* ceiling division */
  d->nslots = ((d->nslots*sizeof(VAL_TYPE) + (GPAGE_SIZE-1)) / GPAGE_SIZE )
            * GPAGE_SIZE;
  d->nslots = d->nslots / sizeof(VAL_TYPE);
    
  d->nelements = 0; /* how many elements currently in the dictionary */
  d->nelements_asked_to_be_inserted = 0;

  d->nprobes_insert=0;
  d->nprobes_lookup=0;

  // the extra d->nslots_per_bucket seems to supress the address sanitizer
  // error, however, i am not sure why since all accesses are < nslots.
  /* d->values = (VAL_TYPE*) aligned_alloc(ALIGNMENT, */
  /* 				      (nslots+d->nslots_per_bucket)*(sizeof(VAL_TYPE))); */
  d->values = (VAL_TYPE*) aligned_alloc(GPAGE_SIZE,
					(d->nslots)*(sizeof(VAL_TYPE))
					+ GPAGE_SIZE);
             /* address sanitizer complains wihtout the above addition */

  madvise(d->values,
	  GPAGE_SIZE,
	  MADV_HUGEPAGE);
  
  /* d->values = (VAL_TYPE*) malloc((nslots)*(sizeof(VAL_TYPE))); */

  // Ensure keys are zeroed
  #pragma omp simd
  for (size_t i = 0; i < nslots; ++i) 
     d->values[i] = 0; // 0 if it is not occupied

  return d;
}

inline void dict_free(dict* d){
  free(d->values);
}
size_t dict_memory(size_t nelements){
  /// return memory estimation of the dictionary size in BYTES
  int nslots_per_bucket = SIMD_LEN; // we store 32 bits per value
  size_t nslots = nelements;
  nslots = nslots + (-nslots % nslots_per_bucket);
  size_t estimate = nslots*(sizeof(VAL_TYPE)) + sizeof(dict);
  
  return estimate;
}







int dict_add_element_to(dict* d, u8* state){
  // =========================================================================+
  // INPUTS:                                                                  |
  // `*d`:  dictionary that will store state as an element                    |
  // `*state`: element to be stored in *d in the form                         |
  // -------------------------------------------------------------------------+

  /* how many bytes do we need to index the buckets */
  const int idx_size =   (int) ceil((log2(NSLOTS_MY_NODE) - log2(d->nslots_per_bucket))
				   /8.0) ;

  ++(d->nelements_asked_to_be_inserted);


  /// Use linear probing to add element to the array d->values
  // get bucket number, recall keys[nbuckets*nslots_per_bucket
  u64 idx = 0;
  memcpy(&idx, state, idx_size);


  /* get the bucket number and scale the index */
  idx = (idx % d->nbuckets) * d->nslots_per_bucket;
  VAL_TYPE val = 0;

  memcpy(&val,
	 &state[idx_size],
	 VAL_SIZE_BYTES );

  /* 0 means empty, we have to ignore zero values  */
  if (val == 0) return 0;

  REG_TYPE dict_slots_simd;
  REG_TYPE value_simd = SIMD_SET1_VALTYPE(val); // (val, val, ..., val) 
  REG_TYPE zero_simd = SIMD_SET1_VALTYPE(0); // (val, val, ..., val) 
  int found_good_slot = 0;

  int zero_idx = 0; /* where to put the element in the dictionary */
  
  // linear probing: examing 16 elements at once using simd
  for (int i=0;
       i<(NPROBES_MAX/SIMD_LEN); /* only Multiples of SIMD_LEN*/
       ++i) {

    
    dict_slots_simd = SIMD_LOAD_SI((REG_TYPE*)  &(d->values[idx]));

    /* is there is a zero */
    found_good_slot = SIMD_CMP_VALTYPE(zero_simd, dict_slots_simd)
                    | SIMD_CMP_VALTYPE(value_simd, dict_slots_simd);

    if (found_good_slot){ /* either it's zero or we have the same value */
      /* the index of the first zero in the simd vector */
      zero_idx = __builtin_ctz(found_good_slot);
      /* put the val in the first zero index  */
      d->values[idx + zero_idx] = val;
      ++(d->nelements);
      return 1;
    }

    idx += d->nslots_per_bucket;
    // reduce mod n->slots //
    if (idx >= d->nslots ) 
      idx = 0;
  }

  return 0; // no element has been added
}








int dict_has_elm(dict *d, u8 *state)
{
  // =========================================================================+
  // returns 1 if state is found in d, 0 otherwise                            |
  // This dictionary is unusual:                                              |
  // -------------------------------------------------------------------------|
  // INPUTS:                                                                  |
  // `*d`:  dictionary that will store state as an element                    |
  // `*state`: element to be looked up  in *d in the form                     |
  // -------------------------------------------------------------------------+

  /* how many bytes do we need to index the buckets */
  const int idx_size =  (int) ceil((log2(NSLOTS_MY_NODE) - log2(d->nslots_per_bucket))
				   /8.0) ;

  u64 idx = 0;
  memcpy(&idx, state, idx_size);

  idx = (idx % d->nbuckets) * d->nslots_per_bucket;

  VAL_TYPE val = 0;
  memcpy(&val,
	 &state[idx_size],
	 VAL_SIZE_BYTES );

  /* 0 means empty, we have to ignore zero values  */
  if (val == 0) return 0;
  

  int is_key_found = 0;
  REG_TYPE dict_keys_simd;// = _mm256_loadu_si256((__m256i*)  &(d->values[h]));
  REG_TYPE lookup_key_simd = SIMD_SET1_VALTYPE(val); // (val, val, ..., val) 
  
  // loop at most NPROBES_MAX/SIMD_LEN since we load SIMD_LEN
  // elements from dictionary each loop.
  for (size_t i=0; i< (int) (NPROBES_MAX/SIMD_LEN); ++i) {
    // we are relying that the val != 0, Pr(val !=0 ) = 1/(2^32)
    // linear probing

    // get new fresh keys from one bucket
    dict_keys_simd = SIMD_LOAD_SI((REG_TYPE*)  &(d->values[idx]));


    // -------------------------------------------------//
    //                   TEST 1:                        //
    /* is one the slots empty? then no point of probing */
    //- ------------------------------------------------//
    // no need for the commented instruction since we only check the first value
    // comp_vect_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);
    // copy the value of the first slot of the bucket, check is it 0?
    if (d->values[idx] == 0)
      return 0;
    
    // -----------------------------------------------//
    //                   TEST 2                       //
    /*  Does key equal one of the slots?              */
    //------------------------------------------------//
    /* 0 no value found, otherwise a value was found  */
    is_key_found = SIMD_CMP_VALTYPE(lookup_key_simd, dict_keys_simd);

    if (is_key_found)
      return 1; /* we will hash the whole message again */


    // update the index for keys load
    idx += d->nslots_per_bucket; // move to the next bucket

    if (idx >= d->nslots)
      idx = 0;

    
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
	   "key = 0x%016x\n",
	   b,
	   d->values[b]);
  }
}
 
