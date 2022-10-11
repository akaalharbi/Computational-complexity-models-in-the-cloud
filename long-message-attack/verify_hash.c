/// The purpose of this file is to check the correctness of the results produced
/// by `long_meessage_attack.c`, namely the messages in the directory
/// `messages/` indeed collide at index `idx` mentioned in the folder
/// statistics. for simplicity we will manually provide `idx` now. We will think
/// later how to automate this.

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "sha256-x86.h"
#include "shared.h"
#include "util/util_char_arrays.h"
int is_there_duplicate = 0;

size_t where_collides(const unsigned char rM[64], int output_size_bits){
  /// Method 1
  /// Given any message rM, hash a long message of zeros M till we find a
  /// collision, then return the shortest length of M that collides with rM.

  size_t idx = 0;


  // for zero messages
  unsigned char M0[64] = {0};
  uint64_t digest_M0[2] = {0}; // store digest here

  uint32_t state0[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };


  uint32_t state_rM[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };

  uint64_t digest_rM[2] = {0};
  truncate_state32bit_get_digest(digest_rM, state_rM, output_size_bits);

  while (1){
    sha256_process_x86_single(state0, M0);
    truncate_state32bit_get_digest(digest_M0, state0, output_size_bits);
    if (digest_M0[0] == digest_rM[0] && digest_M0[1] == digest_rM[1])
      break;
    ++idx;
  }

  return idx;
}


size_t collides_at(const unsigned char rM[64], int output_size_bits, uint64_t idx){
  /// Method 1
  /// Given any message rM, hash a long message of zeros M till we find a
  /// collision, then return the shortest length of M that collides with rM.


  // for zero messages
  unsigned char M0[64] = {0};
  uint64_t digest_M0[2] = {0}; // store digest here

  uint32_t state0[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };


  uint32_t state_rM[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };

  uint64_t digest_rM[2] = {0, 0};
  sha256_process_x86_single(state_rM, rM);
  truncate_state32bit_get_digest(digest_rM, state_rM, output_size_bits);
  
  for (size_t i=0; i<=idx; ++i){
    sha256_process_x86_single(state0, M0);
    truncate_state32bit_get_digest(digest_M0, state0, output_size_bits);
    printf("n=%d, digestrM[0]=%lx, digestM0=%lx, i=%lu\n",output_size_bits, digest_rM[0], digest_M0[0], i);
   

  }
  if (digest_M0[0] == digest_rM[0] && digest_M0[1] == digest_rM[1]){

    return 1;
  }

  return 0;
}


int main(int argc, char* argv[]){
  

  if (argc == 3){
  char* file_name = argv[1];
  FILE* fp = fopen(file_name, "rb");
  int output_size_bits = atoi(argv[2]);
  unsigned char rM[64];
  fread(rM, 64, 1, fp);
  
  
  puts("searching for collision index");  
  size_t idx = where_collides(rM, output_size_bits);

  printf("collides at %lu\n", idx);
  }

  if (argc == 4) {
  char* file_name = argv[1];
  FILE* fp = fopen(file_name, "rb");
  int output_size_bits = atoi(argv[2]);
  size_t at_index = atoi(argv[3]);
  unsigned char rM[64];
  fread(rM, 64, 1, fp);
  fclose(fp);

  puts("message=");
  for (int i = 0; i<64; ++i)
    printf("%2x, ", rM[i]);
  puts("");
  
  FILE* fw = fopen("log/collision_checks", "a"); 

  printf("%d check if it collides at the index %lu\n", output_size_bits, at_index); 
  int collides = collides_at(rM, output_size_bits, at_index);
  fprintf(fp, "%d\n", collides);
  fclose(fw);
  printf("collides?  %d\n", collides);
  size_t idx = where_collides(rM, output_size_bits);
  printf("collides at %lu\n", idx);
  puts("---------");
  }


  else {
    puts("usage:\n1-\n"
	 "./verify_hash message_path output_size_bits");

    puts("2-:\n"
	 "./verify_hash message_path output_size_bits collision_idx");
    
    return 0;
  }


  return 0;
}

