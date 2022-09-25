#ifndef UTIL_CHAR_ARRAYS
#define UTIL_CHAR_ARRAYS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/random.h>

int cmp_arrays(char* array1, char* array2, size_t len);
void print_char(char *l, size_t len);
unsigned char *long_message_zeros(size_t n_of_bits);
unsigned char *create_radom_byte_array(int n_of_bytes);
void fill_radom_byte_array(unsigned char* A, int n_of_bytes);
void truncate_array(unsigned char* A, size_t size_A, size_t total_out_bits);

#endif
