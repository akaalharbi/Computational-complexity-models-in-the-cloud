#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "hash.h"
#include "numbers_shorthands.h"
#include "util_char_arrays.h"
#include "util_files.h"
int check_hashes_interval(const WORD_TYPE state_befoe[NWORDS_STATE],
			  const WORD_TYPE state_after[NWORDS_STATE],
			  const CTR_TYPE ctr){


  u8 M[HASH_INPUT_SIZE] = {0};
  /* We don't want to modify arguments */
  WORD_TYPE state[NWORDS_STATE];
  memcpy(state, state_befoe, HASH_STATE_SIZE);
  
  /* register counter in M */
  ((u64*)M)[0] = ctr;

  for (size_t i=0; i<INTERVAL; ++i) {

    hash_single(state, M);
    ++( ((u64*)M)[0]) ;

    if (memcmp(state, state_after, HASH_STATE_SIZE) == 0){
      printf("i=%lu\n",i );
      print_char((u8*) state, HASH_STATE_SIZE);
      print_char((u8*) state_after, HASH_STATE_SIZE);

    }
  }

  
  return (memcmp(state, state_after, HASH_STATE_SIZE) == 0);
  
}
// @todo discovered a bug in when phase_i rerun again
// it will write the last state twice

int main(int argc, char *argv[])
{

  /* if (argc != 2) */
  /*   return 0; */




  u64 state_number = 3718;

  u64 ctr = state_number*INTERVAL;
  
  WORD_TYPE state_before[NWORDS_STATE];
  WORD_TYPE state_after[NWORDS_STATE];
  
  FILE* fp = fopen("data/states", "r");
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  printf("There are %lu middle states\n", nstates);
  
  fseek(fp, (HASH_STATE_SIZE)*(state_number), SEEK_CUR);

  
  fread(state_before, 1, HASH_STATE_SIZE,  fp);
  fread(state_after, 1,  HASH_STATE_SIZE, fp);

  print_char((u8*) state_before, HASH_STATE_SIZE);
  print_char((u8*) state_after, HASH_STATE_SIZE);
  printf("ctr=%llu, INTERVAL=%llu\n", ctr, INTERVAL);
  puts("Going to check...");

<<<<<<< HEAD
  
  int non_equal_bytes = check_hashes_interval(state_before,
					      state_after,
					      ctr);
  int is_corrupt = 
=======
  int nbytes_non_equal = check_hashes_interval(state_before,
					       state_after,
					       ctr);

  int is_corrupt = (0 != nbytes_non_equal);
>>>>>>> 02a8764b17fd59c1e32605790bdc240c8260083d

  printf("Found a corrupt state? %d\n", is_corrupt);

  fclose(fp);
  return 0;
  
}
