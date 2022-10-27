#include <stdint.h>
#include "intel-headers/arch_avx_type1.h"



typedef struct {
        DECLARE_ALIGNED(uint32_t digest[SHA256_DIGEST_SZ], 32);
        const uint8_t *data_ptr[AVX512_NUM_SHA256_LANES];
} SHA256_ARGS;
