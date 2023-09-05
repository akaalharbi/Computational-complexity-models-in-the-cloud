#include "config.h"



#ifndef LONG_MESSAGE_ATTACK_COMMON
#define LONG_MESSAGE_ATTACK_COMMON
#include <stdio.h>
#include <stdlib.h>

int is_dist_state(u8 state[HASH_STATE_SIZE]);
void print_attack_information();
u32 to_which_server(u8 MState[HASH_STATE_SIZE]);
int is_dist_msg(u8 M[HASH_INPUT_SIZE]);
int is_dist_digest(u8 state[N]);
void transpose_state(u32 dest[restrict 16*8],
		     u32 src[restrict 16*8]);

void untranspose_state(u32 dest[restrict 16*8],
		       u32 src[restrict 16*8]);

void copy_transposed_state(u32 *state, u32 *tr_state, int lane);
void copy_transposed_digest(u8 *digest, u32 *tr_state, int lane);
u64 n_needed_candidates();

static void* my_malloc(size_t size, const char* file, int line) {
    void* ptr = malloc(size);
    if (ptr) {
        printf("Allocated %zu bytes at %p in %s:%d\n", size, ptr, file, line);
    }
    return ptr;
}

static void my_free(void* ptr, const char* file, int line) {
    printf("Freed memory at %p in %s:%d\n", ptr, file, line);
    free(ptr);
}

#define malloc(size) my_malloc(size, __FILE__, __LINE__)
#define free(ptr) my_free(ptr, __FILE__, __LINE__)



#endif
