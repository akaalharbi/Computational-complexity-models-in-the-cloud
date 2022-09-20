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
#include <endian.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "types.h"
#include "util_char_arrays.h"

/*---------------------------------------------------------*/
///                  UTILITY FUNCTIONS                   ///




// int (*functionPtr)(int,int);

// add_element_to_dictionary(dict *dictionary, char *key, size_t value, size_t input_size)

/*---------------------------------------------------------*/


// add_element_to_dictionary(dict *dictionary, char *key, size_t value, size_t input_size)

void nothing(dict *dictionary, dict_key* key, size_t value, size_t key_size){
  // literally do nothing;
}

void long_message_attack(size_t n_of_bits, double l){
//(size_t n_of_bits, size_t n_of_blocks){    
  /*  Mount long message attack on truncated sha256 */
  /// INPUT:
  ///  - n: compression functions output size in bits
  ///  - l: 2^l = how many message blocks we need to store their intermediate values
  /// todo nblocks instead of n_of_blocks!


  /* init */
  // each 512bit block compressed to `n_of_bits` after hash function call
  // we would like to get `n_of_blocks` blocks

  size_t n_of_blocks = (size_t) ceil(pow(2.0, l));
  //  size_t nelements = (size_t) ceil(pow(2.0, n_of_bits));
  int n_of_bytes = (int) ceil( (float)n_of_bits/8);
  BYTE* M = long_message_zeros(n_of_blocks*512);

  /// --------------- printing ----------------------- ///
  printf("INPUT:- %lu number of bits in the compression function\n", n_of_bits);
  printf("      - %lu is the number of blocks\n", n_of_blocks);
  puts("_____________________________________________");
  /// ------------------------------------------------///
  
  // store values in a dictionary
  // dict_new(#elemnets will be stored, element's size)
  // assuming all element have the same size
  dict* d = dict_new(n_of_blocks, n_of_bytes);
  

  // todo start from here tomorrow
  SHA256_CTX ctx;
  // we will zero the excessive zeros, don't edit
  // BYTE buffer[SHA256_BLOCK_SIZE]; 
  
  sha256_init(&ctx);
  
  /* begin the attack */
  // First phase hash an extremely long message 
  sha256_update(&ctx, M, n_of_blocks*64, n_of_bits, d, dict_add_element_to);

  puts("***   long message has been hashed ***");
  dict_print(d, (int) ceil( (float)n_of_bits/8));
  puts("-----------------------");
  
  // Second phase: hash a random message till a collision is found
  int collision_found = 0;
  BYTE* random_message; // 512 bits
  SHA256_CTX ctx2;
  //  STATE intermediate; // union // @remove
  // we will zero the excessive zeros, don't edit
  sha256_init(&ctx2);
  size_t ctr = 0;
  size_t idx;
  
  while (!collision_found) {
    // create a random message of 64 bytes
    random_message = create_radom_byte_array(64);
    /* printf("trials=%lu\n", ctr); */
    /* print_char((char*)ctx2.state, n_of_bytes); */
    /* puts("------"); */
    sha256_init(&ctx2);
    sha256_transform(&ctx2, random_message, n_of_bits);
    // intermediate.state = ctx2.state; @remove
    
    // test for collision and print the results if successful.
    if (dict_has_key(d, (dict_key *) ctx2.state, n_of_bytes)){
      collision_found = 1;
      idx = dict_get_value(d, (dict_key *) &ctx2, n_of_bytes);
      puts("Found a collision with the following details:");
      printf("#random message trials=%lu, index=%lu, M=",ctr, idx);
      dict_key* intermediate = (dict_key *) ctx2.state;
      print_char(intermediate->bytes, n_of_bytes);
      puts("");
      break; // we don't care about the rest of the loop
    }
    // collision_found = 1;

    ++ctr;
  
  }

  /// this part has been commented temporarily, please remove the comment
  /* // Phase verify we have collision */
  /* // message 1 is the long message M  */
  /* // message 2 is random_message || (M after truncation at index idx) */
  /* BYTE* M2 = long_message_zeros( (n_of_blocks - idx)*512); */
  /* memcpy(M2, random_message, 64); // first 64 bytes is the random message */

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
  /* sha256_update(&ctx2, M2, (n_of_blocks-idx)*64, n_of_bits, d, nothing); */


  /* puts("_______________________________________"); */
  /* puts("Hash message 1:"); */
  /* print_char((char*) ctx.state, n_of_bytes); */

  /* puts("Hash message 2:"); */
  /* print_char((char*) ctx2.state, n_of_bytes); */
  /* puts("_______________________________________"); */
  /* // a union trick to read a differen */


  /* //printf("ctr=%lu, idx=%lu, M=",ctr, idx); */
  /* //print_char(intermediate.sate_as_bytes, n_of_bytes); */
 
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
  long_message_attack(n, l);
  

  return 0;
}



