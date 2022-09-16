#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "supherhash.h"
#include "pstdint.h"
#include <stdlib.h>
#include "dict.h"

void print_char(char* l, size_t len){
  for (size_t i = 0; i<len; ++i)
    printf("%hhu, ", l[i]);
  puts("");
}


char* create_random_array(size_t len){
  char* A = (char *)malloc(sizeof(char)*len);
  int d = 0;
  for (size_t i=0; i<len; ++i){

    d = 0 + rand() / (RAND_MAX / (255 - 0 + 1) + 1);
    // printf("rand=%u\n", d);
    A[i] = (char) d;
    // printf("A[%lu]=%hhu\n", i, A[i]);
  }
    
  
  return A;
}


int main(int argc, char* argv[]){
  int n_of_elements = 10;
  int size_of_key = 8;


  // TEST VECTORS FOR THE DICTIONARY //
  // create arrays of elements, each element is an array of bytes
  char** elements = (char**) malloc(sizeof(char*)*n_of_elements);

  // fill each elements with seemingly random bytes
  for (size_t i = 0; i<n_of_elements; ++i) {
    elements[i] = create_random_array(size_of_key);
  }

  
  // print the elements we have
  for (int i = 0; i<n_of_elements; ++i){
    printf("i=%d\n", i);
    print_char(elements[i], size_of_key);
    puts("----------------------------");
  }
  
  // INIT THE DICITONARY //
  // we assume that dict has equal number of bins to the array
  size_t n_of_bins = n_of_elements; // dictionary can hold more, since a bin can hold undetermined number of elements
  dict* d = dict_new(n_of_bins);

  // FILL THE DICTIONARY //
  for (size_t i = 0; i<n_of_elements; ++i) {
    dict_add_element_to(d, elements[i], i, size_of_key);
  }
  puts("all element has been added to the dictionary");

  // LET'S RETRIEVE ALL ELEMENT THAT WE HAVE INSERTED
  size_t value = 0;
  for (size_t i = 0; i<n_of_elements; ++i) {
    value =  dict_get_value(d, elements[i], size_of_key);
    printf("i=%ld, value=%lu, h=%u\n", i, value, (uint32_t)( SuperFastHash(elements[i],  (size_of_key)) % n_of_bins)) ;
  }

  // PRINT STATISTICS OF COLLISIONS IN THE DICTIONARY
  for (size_t i = 0; i<n_of_bins; ++i) {
    printf("i=%ld, number of elements=%lu\n", i, d->bins[i].n_of_elements);
  }

  printf("sizeof(unsigned char)=%ld\n", sizeof(unsigned char));
  
  puts("seems good!");
  return 0;
}

