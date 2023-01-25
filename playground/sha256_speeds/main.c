
#include "sha256_intel_avx/c_sha256_x16_avx512.h"
#include "sha256_intel_avx/c_sha256_oct_avx2.h"
#include "sha256_ni/sha256_ni.h"


#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/random.h>
#include <sys/types.h>
#include "time.h"
#include "timing.h"
#include "vsha256.h"

#include <math.h>


#define NMSGS 16000000LL


void print_bytes(uint8_t* a, size_t length){
  for (size_t i = 0; i<length; ++i) {
    printf("%02x", a[i]);
  }
}


size_t time_sha_ni(){
  uint8_t* msgs = malloc(64*NMSGS);
  
  getrandom(msgs, 64*NMSGS, 1);
  size_t ctr = 0;
  uint32_t* state_ptr;

  double start = wtime();
  for (size_t i = 0; i<NMSGS; ++i) {
    state_ptr = sha256_single(&msgs[i*64]); 
    ctr += ((state_ptr[0] & 0xFF) == 0) ;
  }

  double elapsed = wtime() - start;
  printf("sha_ni elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes \n",
	 elapsed, NMSGS/elapsed, log2(NMSGS/elapsed));


  free(msgs);
  return ctr;

}


typedef struct {
  uint8_t M[16][64];
} msg_avx;


size_t time_sha_avx2(){
  // place holder for the 16 messages
  msg_avx* msgs = malloc(sizeof(msg_avx)*(NMSGS/8));
  // fill the messages
  for (size_t i = 0; i<(NMSGS/8); ++i) {
    getrandom(msgs[i].M , 16*64, 1);
  }
  

  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double start = wtime(); /* timing */

  uint32_t* state_ptr;
  uint32_t state[8];
  
  for (size_t i = 0; i<(NMSGS/8); ++i) {
    state_ptr = sha256_multiple_oct(msgs[i].M);

    for (int i = 0; i<8; ++i) {
      state[i] = state_ptr[i*16];
    }
    ctr += ((state[3] & 0xFF) == 0) ;
  }

  double elapsed = wtime() - start;
  
  printf("sha_avx2_intel elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes\n",
	 elapsed, NMSGS/elapsed,  log2(NMSGS/elapsed));

  free(msgs);
  return ctr;
}


size_t time_sha_avx512(){
  // place holder for the 16 messages
  msg_avx* msgs = malloc(sizeof(msg_avx)*(NMSGS/16));
  // fill the messages
  for (size_t i = 0; i<(NMSGS/16); ++i) {
    getrandom(msgs[i].M , 16*64, 1);
  }
  

  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double start = wtime(); /* timing */

  uint32_t* state_ptr;
  uint32_t state[8];
  
  for (size_t i = 0; i<(NMSGS/16); ++i) {
    state_ptr = sha256_multiple_oct(msgs[i].M);

    for (int i = 0; i<8; ++i) {
      state[i] = state_ptr[i*16];
    }
    ctr += ((state[3] & 0xFF) == 0) ;
  }

  double elapsed = wtime() - start;
  
  printf("sha_avx512_intel elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes\n",
	 elapsed, NMSGS/elapsed,  log2(NMSGS/elapsed));

  free(msgs);
  return ctr;
}





size_t time_vsha(){

  vsha256_setup();

  u32 AVX_ALIGNED h[8][8];
  vsha256_init(h);
  u32 AVX_ALIGNED msg[16][8] = { };
  vsha256_transform(h, msg);

  size_t ctr = 0;
  double start = wtime();
  uint32_t* state_ptr;


  for (size_t i = 0; i<(NMSGS/8); ++i) {
    vsha256_transform(h, msg);
    state_ptr = (uint32_t*) h[0];
    ctr += ((state_ptr[3] & 0xFF) == 0) ;
  }

  double elapsed = wtime() - start;
  
  printf("vsha elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes \n",
	 elapsed, NMSGS/elapsed,  log2(NMSGS/elapsed));
  return ctr;
}





int main(int argc, char* argv[])
{
 /*  OUTPUT: */
 /*  sha_ni elapsed 2.52sec i.e. 6337210.86 hashes/sec = 2^22.595 hashes  */
 /*  sha_avx_intel elapsed 0.69sec i.e. 23277603.04 hashes/sec = 2^24.472 hashes */
 /*  vsha elapsed 3.53sec i.e. 4538645.35 hashes/sec = 2^22.114 hashes  */

  
  time_sha_ni();
  time_vsha();
  time_sha_avx2(); /* intel implementation of sha256 */
  time_sha_avx512(); /* intel avx512 implementation of sha256 */
  // We're giving huge advantages to vhsa, it computes less tasks
  // but it the slowest!

}
