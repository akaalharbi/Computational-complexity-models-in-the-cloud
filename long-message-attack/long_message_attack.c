// Long message attack on sha256
// we make the following modifications:
// 1- The is simply 10*, i.e. no strengthening just add 1 and enough number of
//    zeros to get a multiple of 512-bit input
// 2- The output of a compression function is truncated to n-bit, i.e. in each
//    block of Merkle-Damgard the the arrow to the right is n bit
// 3- store intermediate values

// TODO
// 1- store intermediate values in  sha256.c fil
// 2- doxygene style
// 3- check the code is sound



#include "sha256.h"
#include "dict.h"
#include <bits/types/struct_timeval.h>
#include <endian.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "util/memory.h"
#include <sys/time.h>

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///




// int (*functionPtr)(int,int);
// add_element_to_dictionary(dict *dictionary, char *key, size_t value, size_t input_size)





int is_there_duplicate = 0;

void nothing(dict *dictionary, dict_key* key, size_t value, size_t key_size){
  // literally do nothing;
}

void long_message_attack(size_t n_of_bits, double l){
//(size_t n_of_bits, size_t n_of_blocks){    
  /*  Mount long message attack on truncated sha256 */
  /// INPUT:
  ///  - n: compression functions output size in bits
  ///  - l: 2^l = how many message blocks we need to store their intermediate values
  ///  -fp: a pointer to file where we will store the benchmark information (time,)
  /// todo nblocks instead of n_of_blocks!


  /// PROCESS THE PRARAMETERS
  is_there_duplicate = 0; // global variable to detect cycles
  size_t n_of_blocks = (size_t) ceil(pow(2.0, l));
  //  size_t nelements = (size_t) ceil(pow(2.0, n_of_bits));
  int n_of_bytes = (int) ceil( (float)n_of_bits/8);

  /* FILE * fp; */
  /* fp = fopen("message1", "w"); */
  /* fwrite(M , 1 , n_of_blocks*64 , fp ); */
  /* fclose(fp); */

  

  /* /// --------------- printing ----------------------- /// */
  /* printf("INPUT:- %lu number of bits in the compression function\n", n_of_bits); */
  /* printf("      - %lu is the number of blocks\n", n_of_blocks); */
  /* printf("Let's computer the memeory the program needs\n"); */
  /* printf("- We need %lu bytes for the long message\n", 64*n_of_blocks); */
  /* printf("- We need %lu bytes for the dictionary\n", sizeof(dict)); */
  /* printf("- + %lu bytes for slots \n", 2*n_of_blocks*sizeof(slot)); */
  /* printf("- There are other minor needs we neglect for now\n"); */
  /* printf("n=%lu, l=%f\n", n_of_bits, l); */
  /* printf("n_of_blocks %lu\n", n_of_blocks); */
  /* printf("- sizeof(dict_key)=%lu\n", sizeof(dict_key)); */
  /* printf("- sizeof(size_t)=%lu\n", sizeof(size_t)); */
  /* printf("- sizeof(int)=%lu\n", sizeof(int)); */
  /* printf("- sizeof(slot)=%lu\n", sizeof(slot)); */
  /* puts("_____________________________________________"); */

  // size of  long_message + dict 
  size_t memory_estimate = 64*n_of_blocks + sizeof(dict) + 2*n_of_blocks*sizeof(slot);
  // we need to take into account also the various pointers in this program
  printf("ESTIMATED MEMEORY %lukb\n", memory_estimate/1000);
  
  
  /// ------- INIT ----------------------------------///
  /// timing
  struct timeval begin, end;
  long seconds = 0;
  long microseconds = 0;
  double elapsed = 0;
  gettimeofday(&begin, 0);

  // store values in a dictionary
  // dict_new(#elemnets will be stored, element's size)
  // assuming all element have the same size
  dict* d = dict_new(n_of_blocks, n_of_bytes);
  
  

  /// -------------- PHASE I ------------------------  ///
  /// First phase hash an extremely long message 
  // create a long message
  BYTE* M = long_message_zeros(n_of_blocks*512);
  // INIT SHA256 
  SHA256_CTX ctx;
  sha256_init(&ctx);
  // n_of_bits is the digest length in bits
  sha256_update(&ctx, M, n_of_blocks*64, n_of_bits, d, dict_add_element_to);
  
  /// proposal to save space by evaluating each item of M each time
  // since the message is just 0, we don't need to store all 0
  /*  BYTE* M = long_message_zeros(512); // todo */
  /* for (size_t i=0; i<n_of_blocks; ++i) { */
  /*   sha256_update(&ctx, M, 64, n_of_bits, d, dict_add_element_to); */
  /* } */


  /* puts("--- long message has been hashed ---"); */
  /* dict_print(d, (int) ceil( (float)n_of_bits/8)); */
  /* puts("-----------------------"); */

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  printf("Phase I time: %.3f seconds.\n", elapsed);
  
  /// -------------- PHASE II --------------------  ///
  // Second phase: hash a random message till a collision is found
  /// check if we have a duplicate first
  if (is_there_duplicate)
    // this corresponds to having cycles while hashing the long message
    // i.e. h(mi ... mj) = h(mi ... mk) where k>i
    return; // todo fill this section.

  /// Print memeory usage ///
  long res_mem  = 0;
  long virt_mem = 0;
  puts("========================================");
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("Phase I     : RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("========================================");
  /// --------------- ///
  
  int collision_found = 0;
  BYTE* random_message = (BYTE *) malloc(sizeof(BYTE)*64); // 512 bits
  SHA256_CTX ctx2;
  //  STATE intermediate; // union // @remove
  // we will zero the excessive zeros, don't edit
  sha256_init(&ctx2);
  size_t ctr = 0;
  size_t idx;

  puts("========================================");
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("Phase I done: RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("========================================");
  
  while (!collision_found) {
    // create a random message of 64 bytes
    fill_radom_byte_array(random_message, 64);



    sha256_init(&ctx2);
    sha256_transform(&ctx2, random_message, n_of_bits);

    
    // test for collision and print the results if successful.
    if (dict_has_key(d, (dict_key *) ctx2.state, n_of_bytes)){
      collision_found = 1;
      idx = dict_get_value(d, (dict_key *) &ctx2.state, n_of_bytes);
      puts("Found a collision with the following details:");
      printf("#random message trials=%lu, index=%lu, M=",ctr, idx);
      dict_key* intermediate = (dict_key *) ctx2.state;
      print_char(intermediate->bytes, n_of_bytes);
      puts("");
      break; // we don't care about the rest of the loop
    }
    ++ctr;

  }
  
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  elapsed = seconds + microseconds*1e-6;
  printf("All time: %.3f seconds.\n", elapsed);
  
  puts("- attack has been done");
  puts("========================================");
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("Finished: RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("========================================");


  // todo move this code to another functoin e.g. check collision
  /* /// Saving the two colliding messages & testing they indeed collide */
  /* // Phase verify we have collision */
  /* // message 1 is the long message M */
  /* // message 2 is random_message || (M after truncation after index idx) */
  /* BYTE* M2 = long_message_zeros( (n_of_blocks - idx + 1)*512); */
  /* memcpy(M2, random_message, 64); // first 64 bytes is the random message */
  /* puts("M2="); */
  /* print_char((char*)M2, 64); */
  
  /* // Save message1 to the file message 1 */
  /* FILE * fp; */
  /* fp = fopen("message1", "w"); */
  /* fwrite(M , 1 , n_of_blocks*64 , fp ); */
  /* fclose(fp); */
  /* // Save message2 to the file message 1 */
  /* fp = fopen("message2", "w"); */
  /* fwrite(M2 , 1 , (n_of_blocks - idx)*64 , fp ); */
  /* fclose(fp); */
  
  /* // hash message 1 and message2 and output their digest */
  /* sha256_init(&ctx); */
  /* sha256_update(&ctx, M, n_of_blocks*64, n_of_bits, d, nothing); */
  /* sha256_init(&ctx2); */
  /* sha256_update(&ctx2, M2, (n_of_blocks - idx + 1)*64, n_of_bits, d, nothing); */

  /* puts("_______________________________________"); */
  /* puts("Hash message 1:"); */
  /* print_char((char*) ctx.state, n_of_bytes); */

  /* puts("Hash message 2:"); */
  /* print_char((char*) ctx2.state, n_of_bytes); */
  /* puts("_______________________________________"); */

  /* /\* // indivdually what M2 hashes to  *\/ */
  /* /\* puts("last check M2 will be hashed to something none trivial"); *\/ */
  /* /\* sha256_init(&ctx2); *\/ */
  /* /\* sha256_transform(&ctx2, M2, n_of_bits); *\/ */
  /* /\* print_char((char*) ctx2.state, n_of_bytes); *\/ */
  /* /\* printf("idx=%lu\n", idx); *\/ */


  // free all the used memory 
  free(M);
  free(random_message);
  free(d->slots);
  free(d);
  /* free(ctx); */
  /* free(ctx2); */

  puts("========================================");
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("After free: RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("========================================");

}




int main(int argc, char* argv []){
  // attack(size_t n_of_blocks, size_t n_of_bits)
  
  int n = 0;
  float l = 0;
  
  if (argc != 3){
    puts("USAGE: ./long_message_attac n l");
    puts("n: 0<n<257, the number of bits in the output of the compression function\n"
	 "l:positive double, 2^l is the number of blocks");
    return 0;
  }
  
  // get the input from the user
  n = atoi(argv[1]);
  l = atof(argv[2]);


  // long_message_attack(size_t n_of_blocks, size_t n_of_bits)


  long res_mem  = 0;
  long virt_mem = 0;
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("Before calling the function RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("============================================================");
  long_message_attack(n, l);
  get_memory_usage_kb(&res_mem, &virt_mem);
  printf("After attack finsihed calling the function RES = %lukb, VmSize = %lukb\n", res_mem, virt_mem);
  puts("============================================================");
  

  return 0;
}



