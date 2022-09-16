#ifndef DICT_LINKED_LIST
#define DICT_LINKED_LIST


#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "supherhash.h"
#include "pstdint.h" // maybe no need for that
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include "util_char_arrays.h"

///-----------------------------------------------------///
///                  data structure                     ///
///-----------------------------------------------------///

typedef struct linked_list {
  // we always have an extra empyt element at the end
  char* key; // the key is a list of bytes
  // in our case the value is just an index, todo make it a generic code
  size_t value;
  struct linked_list *next;
} linked_list;


typedef struct {
  // number of element
  size_t n_of_elements; // how many elements in this bin
  linked_list* start;  // first element in the list
  linked_list* last;  // pointer to the last empty element in the linked list
} bin;


typedef struct {
  // number of bins in the dictionary
  // table of bins
  uint32_t  n_of_bins; // how many bins in the dictionary todo remove this serious limtation
  bin* bins; 
} dict;


//-----------------------------------------------------//
//                    functions                        //
//-----------------------------------------------------//


linked_list *linked_list_new();
linked_list* linked_list_add_element_to(linked_list *current,char *key,size_t value, size_t input_size);

bin *bin_new();

dict *dict_new(size_t n_of_bins);
void dict_add_element_to(dict *dictionary,char *key,size_t value,size_t input_size);
int dict_has_key(dict* d, char* key, size_t key_size);
size_t dict_get_value(dict *d, char *key, size_t input_size);
void dict_print(dict *d, size_t key_size);

#endif
