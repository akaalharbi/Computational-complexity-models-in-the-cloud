#include "util_files.h"
#include "numbers_shorthands.h"
#include <stdio.h>

size_t get_file_size(FILE *fp)
{
  /* return file size in bytes */
  u64 old_position = ftell(fp);
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  /* go back to the old position */
  fseek(fp, old_position, SEEK_SET);
  return size;
} 

