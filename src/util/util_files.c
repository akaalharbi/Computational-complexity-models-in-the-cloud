#include "util_files.h"

size_t get_file_size(FILE *fp)
{
  /* return file size in bytes */
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  rewind(fp);
  return size;
} 

