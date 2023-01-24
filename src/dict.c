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

  d->nelements = 0; /* how many elements currently in the dictionary */
  d->nelements_asked_to_be_inserted = 0;

  d->nprobes_insert=0;
  d->nprobes_lookup=0;

  // the extra d->nslots_per_bucket seems to supress the address sanitizer
  // error, however, i am not sure why since all accesses are < nslots.
  d->values = (VAL_TYPE*) aligned_alloc(ALIGNMENT,
				      (nslots+d->nslots_per_bucket)*(sizeof(VAL_TYPE)));
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
  /// return memory estimation of the dictionary size
  int nslots_per_bucket = SIMD_LEN; // we store 32 bits per value
  size_t nslots = nelements;
  nslots = nslots + (-nslots % nslots_per_bucket);
  size_t estimate = nslots*(sizeof(VAL_TYPE)) + sizeof(dict);
  
  return estimate/1000.0;
}




int dict_add_element_to(dict* d, u8* state){
  // =========================================================================+
  // returns 1 if an element has been added, 0 otherwise                      |
  // This dictionary is unusual:                                              |
  // User have a value in the form:                                           |
  // (dist pts, srvr no) || (L bits) || discard || (VAL_SIZE bits) || discard |
  // Dictionary expects user to pass: (L bits) || discard || (VAL_SIZE bits)  |
  // ------------------------------------------------------------------------ |
  // we don't store (dist pts, server) since they are already determined      |
  // The discarded bits between L and VAL_SIZE are due to the fact we move 1  |
  // at least. Those we choose to forget them. The last disarded bits are     |
  // discarded because they will double the dictionary size if we include em  |
  // -------------------------------------------------------------------------|
  // INPUTS:                                                                  |
  // `*d`:  dictionary that will store state as an element                    |
  // `*state`: element to be stored in *d in the form                         |
  //          (L bits) || discard || (VAL_SIZE bits) more precisely:          |
  //          (L_IN_BYTES bytes)  || (VAL_SIZE bits)                          |
  // issues may aris when VAL_SIZE is larger then what is left in the state   |
  //                                                                          |
  // -------------------------------------------------------------------------+
  static int idx_size = MIN(L_IN_BYTES, N-DEFINED_BYTES-VAL_SIZE_BYTES);

    
  ++(d->nelements_asked_to_be_inserted);


  /// Use linear probing to add element to the array d->values
  // get bucket number, recall keys[nbuckets*nslots_per_bucket
  u64 idx = 0;
  memcpy(&idx, state, idx_size);
  /* printf("initially idx=0x%llx, idx_size=%d\n", idx, idx_size); */

  /* get the bucket number and scale the index */
  idx = (idx % d->nbuckets) * d->nslots_per_bucket;

  
  VAL_TYPE val = 0;
  /* @todo address sanitizer detected stack overflow here  */
  /* fread(stream_pt, sizeof(u8), N-DEFINED_BYTES, fp); */
  // L_IN_BYTES=3, N=7, L=20, VAL_SIZE_BYTES=4, state[6]
  // why do we skpid L_IN_BYTES BYTES 
  /* memcpy(&val, &state[L_IN_BYTES], VAL_SIZE_BYTES); */

  memcpy(&val,
	 &state[idx_size],
	 VAL_SIZE_BYTES );

  
  // linear probing 
  for (int i=0;
       i<(NPROBES_MAX/SIMD_LEN)*SIMD_LEN; /* only Multiples of SIMD_LEN*/
       ++i) {

    
    // found an empty slot inside a bucket
     if (d->values[idx] == 0) { // found an empty slot
      d->values[idx] = val;
      ++(d->nelements); /* successfully added an element */
      return 1;
     }


     ++idx;
    // reduce mod n->slots //
    if (idx >= d->nslots ) /* we forgot the equal sign here */
      idx = 0;
  }
  /* printf("missed entry at idx=%llu\n", idx_old); */
  /* print_u16(&(d->values[idx_old]), NPROBES_MAX); */
  return 0; // element has been added
}



int dict_has_elm(dict *d, u8 *state)
{ // returns 1 if state is found in d, 0 otherwise                            |
  // This dictionary is unusual:                                              |
  // User have a value in the form:                                           |
  // (dist pts, srvr no) || (L bits) || discard || (VAL_SIZE bits) || discard |
  // Dictionary expects user to pass: (L bits) || discard || (VAL_SIZE bits)  |
  // ------------------------------------------------------------------------ |
  // we don't store (dist pts, server) since they are already determined      |
  // The discarded bits between L and VAL_SIZE are due to the fact we move 1  |
  // at least. Those we choose to forget them. The last disarded bits are     |
  // discarded because they will double the dictionary size if we include em  |
  // -------------------------------------------------------------------------|
  // INPUTS:                                                                  |
  // `*d`:  dictionary that will store state as an element                    |
  // `*state`: element to be looked up  in *d in the form                     |
  //          (L bits) || discard || (VAL_SIZE bits)                          |
  // -------------------------------------------------------------------------+
  static int idx_size = MIN(L_IN_BYTES, N-DEFINED_BYTES-VAL_SIZE_BYTES);
  u64 idx = 0;
  memcpy(&idx, state, idx_size);

  idx = (idx % d->nbuckets) * d->nslots_per_bucket;

  VAL_TYPE val = 0;
  memcpy(&val,
	 &state[idx_size],
	 VAL_SIZE_BYTES );






  int is_key_found = 0;
  // it's enough to check if the first element is empty
  // in this version we don't need to check if we hit zero or not.
  /* int empty_bucket = 0; //1 - _mm256_testz_si256(comp_vect_simd, comp_vect_simd); */
  // we can remove one of the above variables 
  //+ todo we need to adjust simd instruction according to the type 
  REG_TYPE dict_keys_simd;// = _mm256_loadu_si256((__m256i*)  &(d->values[h]));
  // val:u32 or val:16 depending on the N and L, as 32 or 16 will be stored
  // and the other bits will be stored as index thus the dependency on L.
  REG_TYPE lookup_key_simd = SIMD_SET1_VALTYPE(val); // (val, val, ..., val) 
  //__m256i zero_vect = _mm256_setzero_si256(); // no need for this with buckets
  REG_CMP_TYPE comp_vect_simd;



  
  // loop at most NPROBES_MAX/SIMD_LEN since we load SIMD_LEN
  // elements from dictionary each loop.
  for (size_t i=0; i< (int) (NPROBES_MAX/SIMD_LEN); ++i) {
        
    // we are relying that the val != 0, Pr(val !=0 ) = 1/(2^32)
    // linear probing

    // get new fresh keys from one bucket
    dict_keys_simd = SIMD_LOAD_SI((REG_TYPE*)  &(d->values[idx]));

    // -----------------------------------------------//
    //                   TEST 1                       //
    /*  Does key equal one of the slots?              */
    //------------------------------------------------//

    comp_vect_simd = SIMD_CMP_VALTYPE(lookup_key_simd, dict_keys_simd);
    #ifndef __AVX512F__ /* avx512 is weird */
    is_key_found = (0 == SIMD_TEST(comp_vect_simd, comp_vect_simd));
    #endif
    
    #ifdef __AVX512F__ /* F for Foundation */
    is_key_found = comp_vect_simd;
    #endif

    
    if (is_key_found) {
      return 1; /* we will hash the whole message again */
    }

    // -------------------------------------------------//
    //                   TEST 2: skipped                //
    /* is one the slots empty? then no point of probing */
    //- ------------------------------------------------//
    // no need for the commented instruction since we only check the first value
    // comp_vect_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);
    // copy the value of the first slot of the bucket, check is it 0?
    // with this new version we don't need to check if the first element is 0
    // for two reasons: 1- The loop size is already defined.
    // 2- if we hit zero but kept moving then worst case scenario we wil have a
    // false positive we just need to increase the number of needed hits.
    


    // Linear probing
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
 
