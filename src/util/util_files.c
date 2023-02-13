#include "util_files.h"
#include <stddef.h>
#include <stdio.h>

size_t get_file_size(FILE *fp)
{
  /* return file size in bytes */
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  rewind(fp);
  return size;
}



void merge(FILE* fp_in, FILE* fp_out)
{
  // Write the content of fp_in to fp_out.

  size_t f_size = get_file_size(fp_in);
  /* 10 MB  */
  size_t chunk_size = 10000000LL;
  uint8_t* buffer = (uint8_t*) malloc(chunk_size);

  
  for (size_t j=0;
       j < (f_size / chunk_size) ;
       ++j)
    {
    fread (buffer, 1, chunk_size, fp_in);
    fwrite(buffer, 1, chunk_size, fp_out);
    }
  
  

  /* register the rest of the of elements less than than 1 MB  */
  size_t  bytes_left = f_size - (f_size / chunk_size) * chunk_size;

  fread (buffer, 1, bytes_left, fp_in);
  fwrite(buffer, 1, bytes_left, fp_out);

  free(buffer);
}



// @todo all codes below are untestedx
void write_file_to_another(FILE *fp_in, FILE *fp_out, size_t nbytes)
{
  /* This function doesn't check if there is at least nbytes left in fp_in */
  
  size_t chunk_size = 1000000LL; /* 1KB */
  uint8_t* buffer = (uint8_t*) malloc(chunk_size);
  
  fread (buffer, 1, chunk_size, fp_in);
  fwrite(buffer, 1, chunk_size, fp_out);

  for (size_t j=0;
       j < (nbytes / chunk_size) ;
       ++j)
    {
    fread (buffer, 1, chunk_size, fp_in);
    fwrite(buffer, 1, chunk_size, fp_out);
    }
  
  /* register the rest of the of elements less than than 1 KB  */
  size_t  bytes_left = nbytes - (nbytes / chunk_size) * chunk_size;
  fread (buffer, 1, bytes_left, fp_in);
  fwrite(buffer, 1, bytes_left, fp_out);

  free(buffer);
}


void split_file(FILE *fp, size_t nfiles, size_t word_size)
{
  // Create nfiles each of which has file_size(fp) / word_size words           |
  // i.e.                                                                      |
  // fp = {f_0, f_1, f_2, ..., f_{nfiles-1}}                                   |
  // --------------------------------------------------------------------------+
  
  size_t f_nwords = get_file_size(fp) / word_size;
  size_t nwords_per_file = f_nwords / nfiles;
  FILE* fp_out;
  char file_name[100];

  
  for (size_t f_n = 0; f_n<nfiles; ++f_n) {
    /* create file name string  */
    snprintf(file_name, sizeof(file_name), "data/states/%lu", f_n);
    fp_out = fopen(file_name, "w");

    /* write number of words in file number f_n */
    write_file_to_another(fp, fp_out, nwords_per_file);
    fclose(fp_out);
  }
}
