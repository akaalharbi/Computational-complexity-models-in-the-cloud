#ifndef UTIL_FILES
#define UTIL_FILES

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

size_t get_file_size(FILE *fp);
void merge_file(FILE *fp_in, FILE *fp_out);
void write_file_to_another(FILE *fp_in, FILE *fp_out, size_t nbytes);
void split_file(FILE *fp, size_t nfiles, size_t word_size);

#endif
