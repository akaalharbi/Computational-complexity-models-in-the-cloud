#include "numbers_shorthands.h"
#include "hash.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h> // access functoin

//#include "memory.h" // memory monitor 
//#include <sys/time.h> // moved timing.h 
//#include <assert.h>
#include "config.h"
#include "timing.h"
//#include "types.h" // probably deadweight
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include "util_files.h"




int main(int argc, char* argv[]){
  const size_t hash_state_size = 32;
  size_t old_interval = (1<<30);
  size_t new_interval = (1<<25); 
  size_t enlargen_factor = old_interval/new_interval;
  
  
  FILE* fp = fopen(argv[1], "r"); /* states file is given as an argument */
  size_t fp_size = get_file_size(fp);
  size_t nstates = fp_size / (hash_state_size);
  u8* new_states = malloc(fp_size*enlargen_factor);
  u8* old_states = malloc(fp_size);
  fread(old_states, fp_size, 1, fp);
  fclose(fp);


  /* start regenerating the long message again */
  for (size_t i=0; i<nstates; ++i) {
    u32  state_priv[8];
    u8 M_priv[HASH_INPUT_SIZE] = {0};
    /* register counter in M */
    ((u64*)M_priv)[0] = i*old_interval;

    /* get the middle state we should start from */
    memcpy(state_priv,
	   &old_states[hash_state_size*i],
	   hash_state_size);

    for (size_t j = 0; j<old_interval; j++){
      hash_single(state_priv, M_priv);
      ++( ((u64*)M_priv)[0]) ;

      if (j % new_interval == 0){
	memcpy(state_priv,
	       &new_states[i*enlargen_factor*hash_state_size
			   + j*hash_state_size],
	       hash_state_size);
	
      }
    }

  }


  
  
  return 0;
}
