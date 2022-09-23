/// A collection of useful functions with char arrays
/// extracted while making the long message attack

#include "util_char_arrays.h"
#include <sys/random.h>

int cmp_arrays(char* array1, char* array2, size_t len){

  /// memcmp gives strange results
  for (size_t i = 0; i<len; ++i){
    if (array1[i] != array2[i])
      return 0;
  }
  
  return 1;
}


void print_char(char* l, size_t len){
  printf("0x");
  for (size_t i = 0; i<len; ++i)
    printf("%02x",(unsigned char) l[i]);
  puts("");
}


unsigned char* long_message_zeros(size_t n_of_bits){
  /*     DESCRIPTION        */
  /// Create array of zeros
  /// input: n_of_bits: how many bits the output size should been
  /// output: array of 0 bytes that has `n_of_bits` bits
  /// e.g. input := 9, we can accommedate 9 bits in two bytes
  size_t n_of_bytes = (size_t) ceil(n_of_bits/8.0);
  unsigned char* A = (unsigned char *) malloc(sizeof(unsigned char) * n_of_bytes);
  
  for (size_t i = 0; i<n_of_bytes; ++i)
    A[i] = 0;

  return A;
}


unsigned char* create_radom_byte_array(int n_of_bytes){
  /* Create seemingly a random byte array with total_n_of_bits */
  /// INPUT: how many bytes
  /// OUTPUT: array with ceil(total_n_of_bytes) entries
  ///         the last entry doesn't necessarily use all the 8 bits
  unsigned char* A = (unsigned char *)malloc(sizeof(unsigned char)*n_of_bytes);

  /* int d = 0; */
  /* for (size_t i=0; i<n_of_bytes; ++i){ */
  /*   d = 0 + rand() / (RAND_MAX / (255 - 0 + 1) + 1);    */
  /*   A[i] = (unsigned char) d; */
  /* } */
  int returned_bytes = getrandom(A, n_of_bytes, 1);
  ++returned_bytes; // dummy operation to avoid not used warning

  return A;
}



void fill_radom_byte_array(unsigned char* A, int n_of_bytes){
  /* Create seemingly a random byte array with total_n_of_bits */
  /// INPUT: how many bytes
  /// OUTPUT: array with ceil(total_n_of_bytes) entries
  ///         the last entry doesn't necessarily use all the 8 bits
  

  /* int d = 0; */
  /* for (size_t i=0; i<n_of_bytes; ++i){ */
  /*   d = 0 + rand() / (RAND_MAX / (255 - 0 + 1) + 1);    */
  /*   A[i] = (unsigned char) d; */
  /* } */
  int returned_bytes = getrandom(A, n_of_bytes, 1);
  ++returned_bytes; // dummy operation to avoid not used warning

}
