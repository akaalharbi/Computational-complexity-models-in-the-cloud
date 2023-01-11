#include "config.h"


void print_attack_information();

u32 to_which_server(u8 MState[NWORDS_DIGEST * WORD_SIZE]);
u32 to_which_server(u8 MState[NWORDS_DIGEST * WORD_SIZE]);

void find_hash_distinguished(u8 M[HASH_INPUT_SIZE], /* in, out*/
			     WORD_TYPE Mstate[NWORDS_STATE], /* out*/
			     CTR_TYPE* ctr, /* in, out */
			     const size_t dist_test /* in */);

