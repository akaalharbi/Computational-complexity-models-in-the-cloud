


#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/random.h>
#include <sys/types.h>


#include "sha256_intel_avx/c_sha256_oct_avx2.h"
#include "sha256_ni/sha256_ni.h"


#include <math.h>


#define NMSGS 16000000LL


void print_bytes(uint8_t* a, size_t length){
  for (size_t i = 0; i<length; ++i) {
    printf("%02x", a[i]);
  }
}


void print_digest_1(uint32_t* state){
  uint8_t* state_ptr_u8;
  
  for (int i = 0; i<8; ++i) {
    state_ptr_u8 = (uint8_t*) &state[i*16];
    print_bytes(state_ptr_u8, 4);
  }
  puts("");
}



void print_digest_k(uint32_t* state, int k){
  uint8_t* state_ptr_u8;
  for (int i = 0; i<8; ++i) {
    state_ptr_u8 = (uint8_t*) &state[k + i*16];
    print_bytes(state_ptr_u8, 4);
  }
  puts("");
}


void cpy_transposed_state(uint32_t* tr_state, uint32_t* state, int lane){
  for (int i = 0; i<8; ++i) {
    state[i] = tr_state[lane + i*16];
  }
}



int main(int argc, char* argv[])
{

  uint32_t tmp_state[8] = {0};
  uint8_t* buff = (uint8_t*) malloc(64*16);
  uint8_t msgs[16][64] = {0};
  uint32_t * state = sha256_multiple_8(msgs);
  uint8_t msg[64] = {0};
  uint32_t* state_single;
  int equal = 0;
  int nhashes_dont_agree = 0;
  
for (size_t j = 0; j<100; ++j) {
  getrandom(buff, 64*16, 1);
  memcpy(msgs, buff, 64*16);

  sha256_multiple_8(msgs);


  // go over the 16 messages, hash them and compare them with the lane
  for (int lane=0; lane<8; ++lane) {
    memcpy(msg, &buff[64*lane], 64);
    state_single = sha256_single(msg);
    cpy_transposed_state(state, tmp_state, lane);
    equal = ( 0 == memcmp(tmp_state, state_single, 32));

    if (equal != 1){
      print_bytes(msg, 64);
      puts("at the above message the two hashing methods do not agree!");
      ++nhashes_dont_agree;
    }
  }
 }

 printf("#hashes where two methods of hashing don't agree = %d\n",
	nhashes_dont_agree);

}
