
/// source : https://hpcf.umbc.edu/general-productivity/checking-memory-usage/


#ifndef MEMORY_MONITOR
#define MEMORY_MONITOR
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb);
void print_memory_usage(char* txt);
#endif

