#include "include/dict.h"
#include "include/sha256.h"
#include <stddef.h>



int main(){
  //for (size_t k=25; k<26; k++) {
  size_t k = K_DICT_SIZE; //25;
    size_t dict_size = 1LL<<k;
  
    dict* d = dict_new(dict_size);

    BYTE M[64] = {0}; // long_message_zeros(n_of_blocks*512);

    // store the hash value in this variable
    uint64_t digest[2] = {0, 0};
    uint32_t state[8] = {
      0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    for (size_t i = 0; i<dict_size; i++) {
      sha256_single(state, M);
      truncate_state32bit_get_digest(digest, state, 128);
      dict_add_element_to(d, digest);
    }

    size_t scan_sum[max_seq_length] ={0};
  
    // printf("r_seq%lu = [", k);
    printf("r_seq = [");
    for (int i=0 ; i<max_seq_length; i++) {
      printf("%lu, ",  d->r_sequence[i]);
    }
    puts("]");

    scan_sum[0] = d->r_sequence[0];
    //printf("{%d: %lu}, ", 0, scan_sum[0]);

    for (int i=1 ; i<max_seq_length; i++) {
      scan_sum[i] = scan_sum[i-1] + d->r_sequence[i];
      //printf("{%d: %lu}, ", i, scan_sum[i]);
    }
    
    //printf("mean%lu=%f\n", k,  scan_sum[0] / ( (double) scan_sum[max_seq_length-1])  );
    printf("mean=%f\n",  scan_sum[0] / ( (double) scan_sum[max_seq_length-1])  );

    free(d->keys);
    free(d);
    //}

}
