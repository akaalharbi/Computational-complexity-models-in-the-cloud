#include "config.h"



#ifndef LONG_MESSAGE_ATTACK_COMMON
#define LONG_MESSAGE_ATTACK_COMMON

int is_dist_state(u8 state[HASH_STATE_SIZE]);
void print_attack_information();
u32 to_which_server(u8 MState[HASH_STATE_SIZE]);
int is_dist_msg(u8 M[HASH_INPUT_SIZE]);
void cpy_transposed_state(u32 *Mstate, u32 *tr_state, int lane);

#endif
