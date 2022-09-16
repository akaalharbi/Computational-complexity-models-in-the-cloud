// planning: SERIOUS LIMITATION ONLY 32 BIT DICTIONARY 
//  test the dictionary implementation
// parameters: - size of the dictionary
// hash function: Bernstein hash xor
// methods: - add element: takes a list of bytes and its length
//          - remove element: no need for that now
//


/* #include <stddef.h> */
/* #include <stdint.h> */
/* #include <stdlib.h> */
/* #include "supherhash.h" */
/* #include "pstdint.h" // maybe no need for that */
#include "dict.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// todo refactor this code as c++ code

/* // moved to the header file */
/* //-----------------------------------------------------// */
/* //                  data structure                     // */
/* //-----------------------------------------------------// */

/* typedef struct linked_list { */
/*   // we always have an extra empyt element at the end */
/*   char* key; // the key is a list of bytes */
/*   // in our case the value is just an index, todo make it a generic code */
/*   size_t value; */
/*   struct linked_list *next; */
/* } linked_list; */


/* typedef struct { */
/*   // number of element */
/*   size_t n_of_elements; // how many elements in this bin */
/*   linked_list* start;  // first element in the list */
/*   linked_list* last;  // pointer to the last empty element in the linked list */
/* } bin; */


/* typedef struct { */
/*   // number of bins in the dictionary */
/*   // table of bins */
/*   uint32_t  n_of_bins; // how many bins in the dictionary todo remove this serious limtation */
/*   bin* bins;  */
/* } dict; */


//-----------------------------------------------------//
//                       methods                       //
//-----------------------------------------------------//

linked_list* new_linked_list(){
  linked_list* list = (linked_list *) malloc(sizeof(linked_list));

  //list->key = 0;
  // list->key = (char *) malloc(sizeof(char));
  list->value = 0;
  list->next = 0;
    
  return list;
}



/* bin* new_bin(){ */
/*   bin* slot = (bin *) malloc(sizeof(bin)); */

/*   slot -> n_of_elements = 0; */
/*   slot -> start = new_linked_list(); */
/*   slot -> last = slot->start; */
  
/*   return slot; */
/* } */



dict* dict_new(size_t n_of_bins){
  dict* dictionary = (dict *) malloc(sizeof(dict));
  dictionary->bins =  (bin *) malloc(sizeof(bin) * n_of_bins);
  dictionary->n_of_bins = n_of_bins;

  for (size_t i = 0; i<n_of_bins; ++i) {
    dictionary->bins[i].n_of_elements = 0;
    dictionary->bins[i].start = new_linked_list();
    dictionary->bins[i].last = dictionary->bins[i].start;
    
  }

  
  return dictionary;
}

linked_list* linked_list_add_element_to(linked_list* current, char* key, size_t value, size_t input_size){
  linked_list* next = new_linked_list();

  current->next = next;
  current->key = key;
  //memcpy(current->key, key, input_size);
  current->value = value;
  
  return next;
}

void dict_add_element_to(dict* dictionary, char* key, size_t value, size_t input_size){
  // add (key, value) to the dictionary
  // in this primitive design key is an array of bytes, thus input_size is how many bytes
  // value should be an index of some array thus it has the type size_t
  // todo if we entered a key twice, there will be two entries :(

  // dict = {bin1, bin2, ..., bin_n}
  // find which bin to put it in
  uint32_t h = SuperFastHash(key, input_size);
  // mod is necessary to stay within bounds of dict
  h =(uint32_t) (h % dictionary->n_of_bins); 

  // locate where to place (key, value) within the bin_h
  // dictionary->bins[h];
  bin* slot = &dictionary->bins[h];
  // add (key, value) in the last node of linked list
  slot->last->value = value;

  slot->last->key = (char*) malloc(sizeof(char)*input_size);
  memcpy(slot->last->key, key, input_size);
  slot->n_of_elements += 1;

  // then add a new empty element in the linked list
  slot->last->next = new_linked_list();
  slot->last = slot->last->next;

}



size_t dict_get_value(dict* d, char* key, size_t input_size){
  // we first find to which bin the element should be in
  uint32_t h = SuperFastHash(key, input_size);
  h =(uint32_t) (h % d->n_of_bins);

  
  bin current_slot = d->bins[h];
  linked_list* node = current_slot.start;

  // if the bin is empty todo we should raise an error!
  if (current_slot.n_of_elements == 0)
    return -1;

  // -k1->k2->k3->...->kn
  // we walk through each key ki, then we compare it with input
  char* found_key = node->key;
  for (size_t i=0; i<current_slot.n_of_elements; ++i){ 
    if (cmp_arrays(found_key, key, input_size)) 
      return node->value; // we found key

    node = node->next; // ki->k{i+1}
    found_key = node->key; 
  }
  
  return -1; // we should throw an error here, todo
}


void dict_print(dict* d, size_t key_size){
  size_t n_of_bins = d->n_of_bins;

  for (size_t b=0; b<n_of_bins; ++b) {
    bin* slot = &d->bins[b];
    linked_list* node = slot->start;

    printf("bin=%lu, n_of_elements=%lu\n", b, slot->n_of_elements);
    for (size_t j=0; j<slot->n_of_elements; ++j) {
      printf("(val=%lu, key={0x", node->value);// todo print the key as well!
      // printing the key
      for (size_t k=0; k<key_size; ++k)
	printf("%x",(unsigned char) node->key[k]);
      printf("}), ");
      node = node->next;
    }
    puts("\n______________________________");
  }
  
}


int dict_has_key(dict* d, char* key, size_t key_size){
  /// 1 if the dictionary has `key`
  /// 0 otherwise

  
  // we first find to which bin the element should be in
  uint32_t h = SuperFastHash(key, key_size);
  h =(uint32_t) (h % d->n_of_bins);

  
  bin current_slot = d->bins[h];
  linked_list* node = current_slot.start;

  // if the bin is empty todo we should raise an error!
  if (current_slot.n_of_elements == 0)
    return 0;

  // -k1->k2->k3->...->kn
  // we walk through each key ki, then we compare it with input
  char* found_key = node->key;
  for (size_t i=0; i<current_slot.n_of_elements; ++i){ 
    if (cmp_arrays(found_key, key, key_size)) 
      return 1;

    node = node->next; // ki->k{i+1}
    found_key = node->key; 
  }
  
  return 0; // we should throw an error here, todo
 
}

