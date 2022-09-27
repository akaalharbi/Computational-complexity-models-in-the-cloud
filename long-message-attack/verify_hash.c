/// The purpose of this file is to check the correctness of the results produced
/// by `long_meessage_attack.c`, namely the messages in the directory
/// `messages/` indeed collide at index `idx` mentioned in the folder
/// statistics. for simplicity we will manually provide `idx` now. We will think
/// later how to automate this.

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "sha256.h"
#include "shared.h"
int is_there_duplicate = 0;

size_t where_collides(const unsigned char rM[64], int output_size_bits){
  /// Method 1
  /// Given any message rM, hash a long message of zeros M till we find a
  /// collision, then return the shortest length of M that collides with rM.

  size_t idx = 0;
  int does_it_collide = 0;

  // for zero messages
  unsigned char M0[64] = {0};
  uint64_t digest_M0[2] = {0}; // store digest here
  SHA256_CTX ctx;
  sha256_init(&ctx);

  
  SHA256_CTX ctx2;
  sha256_init(&ctx2);
  sha256_transform(&ctx2, rM, output_size_bits);
  uint64_t digest_rM[2] = {0};
  truncate_state_get_digest(digest_rM, &ctx2, output_size_bits);

  while (1){
    sha256_transform(&ctx, M0, output_size_bits);
    truncate_state_get_digest(digest_M0, &ctx, output_size_bits);
    if (digest_M0[0] == digest_rM[0] && digest_M0[1] == digest_rM[1])
      break;
    ++idx;
  }

  return idx;
}



int main(int argc, char* argv[]){
  

  if (argc  != 3){
    puts("usage:\n"
	 "./verify_hash message_path output_size_bits");
    return 0;
  }

  char* file_name = argv[1];
  FILE* fp = fopen(file_name, "rb");
  int output_size_bits = atoi(argv[2]);
  unsigned char rM[64];
  fread(rM, 64, 1, fp);

  puts("searching for collision index");  
  size_t idx = where_collides(rM, output_size_bits);

  printf("collides at %lu\n", idx);
  
  return 0;
}

