#include <stdio.h>
#include <stdint.h>

#ifdef __SIZEOF_INT128__
typedef unsigned __int128      u128;
#endif


void increment_by_one(unsigned __int128* u){
  (*u)++;
}
int main(int argc, char* argv[]){
  /* unsigned char bytes[8] = {255, 2, 3, 4, 250, 5, 6, 7};//{1, 2, 3, 4, 5, 6, 7, 8}; */

  /* uint32_t* word = (uint32_t*) bytes; */
  /* uint64_t* dword = (uint64_t*) bytes; */
  /* // it uses little endian */
  


  /* uint32_t state[8] = { */
  /*   0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, */
  /*   0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 */
  /* }; */

  /* word = (uint32_t*) state; */
  /* dword = (uint64_t*) state; */
  /* printf("As uint32=0x%X, sizeof(bytes)=%lu\n", word[0], sizeof(bytes)); */
  /* printf("As uint32=0x%lX, sizeof(bytes)=%lu\n", dword[0], sizeof(bytes)); */
  /* for (int i=0; i<8; i++) { */
  /*   printf("%02x", ((unsigned char*) state)[i]); */
  /* } */
  /* puts(""); */

  uint64_t a  = (7LL <<32) | 1;
  uint32_t al = (uint32_t) a;
  uint32_t ah = (uint32_t) (a>>32);
  printf("a=%lx, ah=%x, al=%x\n", a, ah, al);

  unsigned __int128 big = -1;
  al = big;
  ah = big >> 64;
  printf(" ah=%x, al=%x\n", ah, al);

  //increment_by_one(&big);
  
  al = big;
  ah = big >> 64;
  printf(" ah=%x, al=%x\n", ah, al);

 uint32_t M[8] = { 0x6a09e667, 0xbb67ae85,
			    0x3c6ef372, 0xa54ff53a,
			    0x510e527f, 0x9b05688c,
			    0x1f83d9ab, 0x5be0cd19 };

 printf("M[0]=%x, M[1]=%x\n", M[0], M[1]);
 u128* offset = (u128*) M;

  offset[0] = big; /* increase the first 128bits by one */

 
 for (int i = 0; i<8; i++) {
   printf("M[%d]=%x, ", i, M[i]); 
 }


 
  return 0;
}
