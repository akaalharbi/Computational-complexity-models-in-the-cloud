#include <stdio.h>
#include <stdint.h>

extern void sha256_block_avx(const void *, void *);


int main(int argc, char* argv[]){
  unsigned char data[64] = {0};
  unsigned char state[32] = {0};

  sha256_block_avx(data, state);

  for (int i=0; i<32; i++) {
    printf("%2x", state[i]);
  }
  puts("");
}
