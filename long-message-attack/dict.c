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
  /// dict v3: we don't store values, we only store 64bit of the key
  ///         in d->keys[k]. When d->keys[k]=0, it means that we have
  ///         we have an empty slot


  dict* d = (dict *) aligned_alloc(ALIGNMENT, (sizeof(dict)));

  
  // Use a defined filling rate type.h, empiricall data shows 0.77
  // is the best
  size_t nslots = (size_t)  ceil((1/FILLING_RATE)*nelements);
  
  // the number of slots is multiple of alignment
  // this is for SIMD. Don't confuse it with sha256 padding
  // multiple of ALIGNMENT BYTES. @update: not sure if this is required!
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
				
  /* d->memory_estimates = (padding_alignment + nslots)*(sizeof(uint64_t)) // keys */
  /*                     + sizeof(size_t)*2 */
  /*                     + sizeof(dict); */
  

  // Ensure the keys are next to each other
  #pragma omp simd
  for (size_t i = 0; i < nslots; ++i) {
     d->keys[i] = 0; // 0 if it is not occupied
  }  
   // printf("- Dictionary of size 0x%lx\n has been initialized\n", nslots);
  return d;
}

inline void dict_free(dict* d){
  free(d->keys);
}


int dict_memory(size_t nelements){
  /// return memory estimation of the dictionary size
  int estimate = (size_t)  ceil((1/FILLING_RATE)*nelements);
  int padding_alignment = (ALIGNMENT/sizeof(uint64_t)) - nelements%( (ALIGNMENT/sizeof(uint64_t)) );
  // how many slots in the dictionary
  estimate = estimate + padding_alignment;

  estimate = estimate*(sizeof(uint64_t)) /* keys */
           + sizeof(dict);

  return estimate/1000.0;
}


void dict_add_element_to(dict* d, uint64_t key[2]){
  /// Use linear probing to add element to the array d->keys
  uint64_t h = key[0]; // for now we assume indices don't need more than 64 bits 
  h = h % d->nslots; // assumption nslots = 2^m; 

  // locate where to place (key, value)
  while (d->keys[h]){ /* ==0; means empty slot */
    // linear probing
    
    ++h;
    if (h >= d->nslots)
      h = 0; //h - d->nslots;
    
    #ifdef NPROBES_COUNT
    ++(d->nprobes_insert);
    #endif
  }

  /// Found an empty slot, update its 
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






//size_t dict_get_values_simd(dict* d, __m256i keys){
void dict_get_values_simd(dict* d, uint64_t keys[4], uint64_t found_keys[4]){
  /// input dict* d,
  /// keys = {k0, k1, k2, ..., kl} 
  /// l is an argument depends on the largest vector lenght available in simd
  /// also, alignment has to be adjusted according to l, in our case 64*l
  ///
  /// output:
  /// values = {v0, v1, v2, ..., vl}
  /// vi = 0 if keys doesn't exist in the dictionary, otherwise return ki!
  /// @todo write about dictionary
  /// comments that have the prefix //+ indicates add programming lines
  /// that correspond to this psuedo-code 
  
  // -------------------- VARIBLES SETUP ------------------------ //
  //+ set component i: 64 bit to be the key passed from input
  __m256i lookup_keys_simd = _mm256_setr_epi64x(keys[0],
						keys[1],
						keys[2],
						keys[3]);
  //+ indices to load from the memory:
  uint64_t indices[4] __attribute__((aligned(32))) = {
      keys[0] % d->nslots,
      keys[1] % d->nslots,
      keys[2] % d->nslots,
      keys[3] % d->nslots
    } ;
  __m256i indices_simd = _mm256_load_si256((__m256i*) indices);
  
  // Load keys from dictionary according to indices: (is gather better?)
  // they will be compared against lookup_keys_simd
  __m256i dict_keys_simd = _mm256_set_epi64x(d->keys[indices[0]],
					     d->keys[indices[1]],
					     d->keys[indices[2]],
					     d->keys[indices[3]]);

  __m256i zero_vect = _mm256_setzero_si256();
  __m256i comp_keys_simd  = _mm256_cmpeq_epi64(lookup_keys_simd, dict_keys_simd);
  __m256i comp_empty_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);

  //+ if we found key i, or we hit an empty slot in the linear probing then there is
  // no point to search further inside the dictionary. In this case, it is same as
  // saying we move 0 step. Thus, component i: 32 bit of steps == 1 if no key nor
  // empty were encountered, 0 otherwise.

  __m256i ones = _mm256_set1_epi64x(1); /* (1, 1, 1, 1) */

  // Currently: 0 if no key were found nor an empty hole,
  // i.e. probe the next element. The 0 will become 1, just be patient
  __m256i steps = _mm256_or_si256(comp_keys_simd, comp_empty_simd);
  // bi = 0 if a key is found or if an empty entry is found
  // bi = 1 if a key were not found nor an empty slot
  // Note:
  // First converts 0 to 1 and vice versa, then AND it with 1
  steps = _mm256_andnot_si256(steps, ones);  ;/* (b0, b1, b2, b3) */  
  indices_simd = _mm256_add_epi64(indices_simd, steps);
  _mm256_store_si256((__m256i*) indices, indices_simd);
  // Note: we can use steps as a stopping indicator
  // If all values are zeros, then it will return 1
  int should_stop = _mm256_testz_si256(steps, steps);

  // @todo modular arithmetic Barret's reduction maybe
  // @todo load from memory using indices

  // @todo extract keys for return
  //+ rewirte the codition to suit the vectorized version
  // @TODO start from here



#ifdef VERBOSE_LEVEL
  #if VERBOSE_LEVEL == 2
  puts("-----------------------");
  puts("Inside dict_get_value ");
  for (int i=0; i<4; i++) 
    printf("k%d=%016lx\n", i, keys[i]);

  for (int i=0; i<4; i++) 
    printf("idx%d=%ld\n", i, indices[i]);
  
  printf("nslots=%lu\n",  d->nslots);
  print_m25i(lookup_keys_simd, "lookup_key_simd");
  print_m25i(comp_keys_simd, "comp_keys_simd");
  print_m25i(dict_keys_simd, "dict_keys_simd");
  print_m25i(comp_empty_simd, "comp_empty_simd");
  print_m25i(steps, "steps");
  printf("should stop=%d\n", should_stop);
  puts("ENTERING WHILE LOOP");
  #endif
#endif
  while (!should_stop){// occupied slo,tmul
    // linear probing
    // get new fresh keys
    //+ change load to be consistent with idx vector

    // update indices according to steps @todo this is ugly
    // for loop does not work here because extract requires const int

    /* indices[0] += _mm256_extract_epi64(steps,  0); */
    /* indices[1] += _mm256_extract_epi64(steps,  1); */
    /* indices[2] += _mm256_extract_epi64(steps,  2); */
    /* indices[3] += _mm256_extract_epi64(steps,  3); */
    // if an index passes the end of the array, go back to the beginning
    #pragma omp simd
    for (int i = 4; i<4; ++i) {
      if (indices[i] >= d->nslots) {
        indices[i] = 0; // reduction mod d->nslots, since the largest step is 1
      }
    }

 
    // Load fresh new keys // use load instead
    dict_keys_simd = _mm256_set_epi64x(d->keys[indices[0]],
				       d->keys[indices[1]],
				       d->keys[indices[2]],
				       d->keys[indices[3]]);


      /// -----------------------------------------------///
     ///                   TEST 1                       ///
    ///------------------------------------------------///
    /* for some i, does ki == dict_keys_simd[i] */
    comp_keys_simd  = _mm256_cmpeq_epi64(lookup_keys_simd, dict_keys_simd);

    ///----------------------------------------------------///
    ///                   TEST 2                           ///
    ///----------------------------------------------------///
    /// is one of the slots empty? then no point of probing for the respective key
    comp_empty_simd = _mm256_cmpeq_epi64(dict_keys_simd, zero_vect);

    ///----------------------------------------------------///
    ///                UPDATE STEPS                        ///
    ///----------------------------------------------------///
    ones = _mm256_set1_epi64x(1); /* (1, 1, 1, 1) */
    // 0xfffffffff if no key were found nor an empty hole, i.e. probe the next element
    steps = _mm256_or_si256(comp_keys_simd, comp_empty_simd);

    // bi = 0 if a key is found or if an empty entry is found
    // bi = 1 if a key were not found nor an empty slot
    // First converts 0 to 1 and vice versa, then AND it with 1
    steps = _mm256_andnot_si256(steps, ones);  ;/* (b0, b1, b2, b3) */  
    indices_simd = _mm256_add_epi64(indices_simd, steps);
    _mm256_store_si256((__m256i*)indices, indices_simd);  
    // Note: we can use steps as a stopping indicator
    // If all values are zeros, then it will return 1
    should_stop = _mm256_testz_si256(steps, steps);

    

#ifdef VERBOSE_LEVEL
    #if VERBOSE_LEVEL == 2
    printf("nslots=%lu\n",  d->nslots);
    print_m25i(lookup_keys_simd, "lookup_key_simd");
    print_m25i(comp_keys_simd, "comp_keys_simd");
    print_m25i(dict_keys_simd, "dict_keys_simd");
    print_m25i(comp_empty_simd, "comp_empty_simd");
    print_m25i(steps, "steps");
    printf("should stop=%d\n", should_stop);
    for (int i=0; i<4; i++) 
      printf("idx%d=%ld\n", i, indices[i]);
    print_m25i(indices_simd, "indices_simd");
    puts("");
    #endif
#endif
    
    
    #ifdef NPROBES_COUNT
    ++(d->nprobes_lookup);
    #endif
  }

  // update found_keys
  _mm256_storeu_si256((__m256i*)indices, dict_keys_simd);
  /* found_keys[0] = _mm256_extract_epi64(dict_keys_simd,  0); */
  /* found_keys[1] = _mm256_extract_epi64(dict_keys_simd,  1); */
  /* found_keys[2] = _mm256_extract_epi64(dict_keys_simd,  2); */
  /* found_keys[3] = _mm256_extract_epi64(dict_keys_simd,  3); */

}


void dict_print(dict* d){
 
  for (size_t b=0; b<(d->nslots); ++b) {
    printf("slot=%lu, "
	   "key = 0x%016lx\n",
	   b,
	   d->keys[b]);
  }
}
 
