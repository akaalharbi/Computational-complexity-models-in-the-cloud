#ifndef UTIL_CHAR_ARRAYS
#define UTIL_CHAR_ARRAYS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int cmp_arrays(char* array1, char* array2, size_t len);
void print_char(char *l, size_t len);
unsigned char *long_message_zeros(size_t n_of_bits);
unsigned char *create_radom_byte_array(int n_of_bytes);

#endif
