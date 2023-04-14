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
  
  
  FILE* fp = fopen(argv[1], "r");
  size_t fp_size = get_file_size(fp);
  size_t nstates = fp_size / (hash_state_size);
  u8* new_states = malloc(fp_size*enlargen_factor);
  fclose(fp);

  for (size_t i=0; i<nstates; ++i) {
    u32  state_priv[8];
    FILE* fp_priv = fopen(argv[1], "r"); /*private to a thread */
    size_t ctr_priv = i*old_interval;

    fseek(fp_priv, (hash_state_size)*(i), SEEK_CUR);
    fread(state_priv, 1, hash_state_size,  fp);
    
    
    
  }
  
  
  return 0;
}
