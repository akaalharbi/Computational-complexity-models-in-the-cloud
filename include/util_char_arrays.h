#ifndef UTIL_CHAR_ARRAYS
#define UTIL_CHAR_ARRAYS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/random.h>
#include <stdint.h>

uint64_t linear_search(uint8_t *key,
		      uint8_t *array,
		      size_t array_len,
                      size_t key_len);

void *linear_search_ptr(uint8_t *key,
			uint8_t *array,
			size_t array_len,
                        size_t key_len);


void print_byte_array(uint8_t* array, size_t nbytes);
void print_u16(uint16_t* l, size_t len);
int cmp_arrays(char* array1, char* array2, size_t len);
void print_char(unsigned char *l, size_t len);
void print_byte_txt(char* txt, unsigned char* a, size_t len);
unsigned char *long_message_zeros(size_t n_of_bits);
unsigned char *create_radom_byte_array(int n_of_bytes);
void fill_radom_byte_array(unsigned char *A, int n_of_bytes,
                           unsigned int *seed);
void fill_radom_byte_array_get_random(unsigned char* A, int n_of_bytes);
void truncate_array(unsigned char* A, size_t size_A, size_t total_out_bits);
void truncate_state32bit_get_digest(uint64_t* dst, uint32_t state[8], int n_of_bits);
#endif
