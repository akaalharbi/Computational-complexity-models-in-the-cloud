// dummy file that will be removed later
#include "numbers_shorthands.h"
#include "hash.h"

#include "dict.h"
#include <bits/types/struct_timeval.h>

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>
#include "config.h"
#include "timing.h"
#include "types.h"
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include <mpi.h>

// MPI sending and receiving tags
#define TAG_DICT_SND 0 /* send digest tag */
#define TAG_RANDOM_MESSAGE 1
#define TAG_SND_DGST 2
#define TAG_MESSAGES_CANDIDATES 3



















// when a process born:

//+ i am a sender:
//+ generate hashes and
//+ keep going

//+ i am a receiver:
//+ load hashes from the file
//+ check how many hashes did i generate?
//+ generate more
//+ how about a deadlock? when a process generates enough hash


//+ i am process with rank NSERVERS (master)
//+ check how many candidates do i have.
//+ enough candidate?
//++ python script try to kill other servers
//+ not enough:
//+ keep listening to other receivers, and record any candidates
//+ update the number of candidate messages 


