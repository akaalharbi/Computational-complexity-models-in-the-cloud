// Long message attack
#include "numbers_shorthands.h"
#include "hash.h"

#include "dict.h"

// deadweight
// #include <bits/types/struct_timeval.h>

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h> // access functoin

//#include "memory.h" // memory monitor 
//#include <sys/time.h> // moved timing.h 
//#include <assert.h>
#include "config.h"
#include "timing.h"
//#include "types.h" // probably deadweight
#include "util_char_arrays.h"
#include "shared.h" // shared variables for duplicate 
#include "memory.h"
#include "util_files.h"
#include "common.h"
#include "dict.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include "c_sha256_avx.h"


//---------------------------- UTILITY FUNCTIONS -------------------------------





// -----------------------------------------------------------------------------
// local function
static void load_file_to_dict(dict *d, FILE *fp)
{
  // ==========================================================================+
  // Summary: Load hashes from the file *fp, and try to store them in dict *d  |
  // Note: dict* d has the right to reject inserting element.                   |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d`: Dictionary that will keep elements from *fp.                        |
  // `*fp` : File contain number of hashes larger than nelements               |
  // ==========================================================================+

  /* Check that file exists, the file comes from external resources */  
  if (!fp){
    puts("I have been given a file that lives in nowhere");
    return; // raise an error instead
  }

  
  size_t nchunks = 10000; 
  size_t ndigests = get_file_size(fp) / N;
  size_t nmemb = (ndigests/nchunks >=  1) ? ndigests/nchunks : 1;

  /* read digests from file to this buffer  */
  u8* digests = (u8*) malloc( N * nmemb * sizeof(u8));


  /*  load one chunk each time */
  for (size_t i = 0; i<nchunks; ++i) {
    fread(digests, N*nmemb, 1, fp);

    /* add them to dictionary */
    for (size_t j=0; j<nmemb; ++j) 
      dict_add_element_to(d,
			  &digests[j*N /* need to skip defined bytes */
				   + DEFINED_BYTES] );
  }

  /* ndigests = nchunks*nmemb + remainder  */
  // read the remainder digests, since ndigests may not mutlipe of nchunks
  u8 stream_pt[N];

  /* add as many hashes as possible */
  while ( !feof(fp) ){
    // use fread with a larger buffer @todo
    fread(stream_pt, sizeof(u8), N, fp);
    /* it adds the hash iff nprobes <= NPROBES_MAX */
    dict_add_element_to(d, &stream_pt[DEFINED_BYTES]);
  }

  free(digests);
  return;
}



typedef struct {
  uint8_t M[16][64];
} msg_avx;


#ifdef  __AVX512F__
size_t time_sha_avx512(){
  // place holder for the 16 messages
  
  #pragma omp parallel
  {
  size_t nmsgs = 1600000LL;
  msg_avx msg;
  getrandom(msg.M , 16*HASH_INPUT_SIZE, 1);
  // fill the messages



  
  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double elapsed = 0;
  double start = wtime(); /* timing */

  for (size_t i = 0; i<(nmsgs/16); ++i) {
    sha256_multiple_x16(msg.M);
  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("thd%d sha_avx512_16way  elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 omp_get_thread_num(),
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);

  }
  return 0;
}



size_t time_sha_one_avx512_other_sha_ni(){
  // place holder for the 16 messages

  #pragma omp parallel
  {
  if (omp_get_thread_num() == 0) {
  size_t nmsgs = 16000000LL;
  msg_avx msg;
  getrandom(msg.M , 16*HASH_INPUT_SIZE, 1);
  // fill the messages



  
  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double elapsed = 0;
  double start = wtime(); /* timing */



  uint32_t* state_ptr;
  
  for (size_t i = 0; i<(nmsgs/16); ++i) {
    state_ptr = sha256_multiple_x16(msg.M);
  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("thd%d sha_avx512_16way  elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 omp_get_thread_num(),
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);

  } else {
    size_t nmsgs = 1600000LL;
    u8 M[64] = {0};
    getrandom(M, 64, 1);
    uint32_t state[8] = {HASH_INIT_STATE};
    double elapsed = 0;
    size_t dummy = 0;
    double start = wtime(); /* timing */
    for (size_t i = 0; i<nmsgs; ++i) {
      hash_single(state, M);
    }

    elapsed = wtime() - start;
    double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);

    printf("thd%d sha_ni elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	   omp_get_thread_num(),
	   elapsed, (nmsgs)/elapsed,
	   log2((nmsgs)/elapsed),
	   hashed_MB/elapsed);
    printf("%lu\n", dummy);
    
    
  }

  }
  return 0;
}


size_t time_sha_avx512_single_thread(){
  // place holder for the 16 messages



  size_t nmsgs = 1600000LL;
  msg_avx msg;
  getrandom(msg.M , 16*HASH_INPUT_SIZE, 1);
  // fill the messages



  
  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double elapsed = 0;
  double start = wtime(); /* timing */


  uint32_t state[8] = {HASH_STATE_SIZE};
  uint32_t* state_ptr = state;
  
  for (size_t i = 0; i<(nmsgs/16); ++i) {
    state_ptr = sha256_multiple_x16(msg.M);
    ctr += ((state[3] & 0xFF) == 0) ;
  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("single sha_avx512_16way  elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);


  return 0;
}
#endif

size_t time_sha_avx256(){
  // place holder for the 16 messages

  #pragma omp parallel
  {
  size_t nmsgs = 8000000LL;
  msg_avx msg;
  getrandom(msg.M , 16*HASH_INPUT_SIZE, 1);
  // fill the messages



  
  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double elapsed = 0;
  double start = wtime(); /* timing */


  uint32_t state[8] = {HASH_INIT_STATE};
  uint32_t* state_ptr = state;
  
  for (size_t i = 0; i<(nmsgs/8); ++i) {
    state_ptr = sha256_multiple_oct(msg.M);
    ctr += ((state[3] & 0xFF) == 0) ;
  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("thd%d sha_8way  elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 omp_get_thread_num(),
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);

  }
  return 0;
}


size_t time_sha_avx256_single(){
  // place holder for the 16 messages



  size_t nmsgs = 8000000LL;
  msg_avx msg;
  getrandom(msg.M , 16*HASH_INPUT_SIZE, 1);
  // fill the messages



  
  size_t ctr = 0; /* dummy variable so that the compiler won't optimize */
  double elapsed = 0;
  double start = wtime(); /* timing */


  uint32_t state[8] = {HASH_INIT_STATE};
  uint32_t* state_ptr = state;
  
  for (size_t i = 0; i<(nmsgs/8); ++i) {
    state_ptr = sha256_multiple_oct(msg.M);
    ctr += ((state[3] & 0xFF) == 0) ;
  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("single sha_8way  elapsed %0.2fsec i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);


  return 0;
}



void bench_long_message_gen()
{
  size_t nmsgs = (1LL<<25);
  
  u8 Mavx[16][HASH_INPUT_SIZE] = {0};
  u32 tr_states[16*8] = {0}; /* same as current_states but transposed */

  getrandom(tr_states, 16*8, 1);
  getrandom(Mavx , 16*HASH_INPUT_SIZE, 1);

  double elapsed = 0;
  double start = wtime(); /* timing */


  
  for (size_t i = 0; i<(nmsgs/NHASH_LANES); ++i) {
      /* hash 16 messages and copy it to tr_states  */    
      #ifdef  __AVX512F__
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states, 0),
	     16*HASH_STATE_SIZE);
      #endif

      #ifndef  __AVX512F__
      #ifdef    __AVX2__
      /* sha256_multiple_oct_tr(Mavx, tr_states); */

      memcpy(tr_states,
	     sha256_multiple_oct_tr(Mavx, tr_states),
	     8*HASH_STATE_SIZE);
      #endif
      #endif

    
    /* update message counters */
    for (int lane = 0; lane<16; ++lane)
      ((u64*) Mavx[lane])[0] += 1;

  }

  elapsed = wtime() - start;
  double hashed_MB = (nmsgs*HASH_INPUT_SIZE) / ((double) 1000000);
  
  printf("regenarate long_message  elapsed %0.2fsec"
	 "i.e. %0.2f hashes/sec = 2^%0.3f hashes, %0.4f MB\n",
	 elapsed, (nmsgs)/elapsed,
	 log2((nmsgs)/elapsed),
	 hashed_MB/elapsed);

 
}


int main(int argc, char *argv[])
{

  print_attack_information();

  printf("Estimated memory per dictionary=%lu Bytes ≈ 2^%0.2f Bytes\n",
	 dict_memory(NSLOTS_MY_NODE),
	 log2(dict_memory(NSLOTS_MY_NODE)));

  #ifdef  __AVX512F__
  /* time_sha_avx512(); */
  #endif
  /* puts("==============================================\n"); */
  #ifdef  __AVX512F__ 
  time_sha_avx512_single_thread();
  #endif
  /* puts("==============================================\n"); */
  /* time_sha_one_avx512_other_sha_ni(); */
  /* puts("==============================================\n"); */
  /* time_sha_avx256(); */
  /* puts("==============================================\n"); */
  /* time_sha_avx256_single(); */
  /* puts("==============================================\n"); */
  bench_long_message_gen();
  puts("==============================================\n");


  return 0;
  printf("sizeof(dict)=%lu bytes\n", sizeof(dict));

  /* Benchmark dictionary query  */
  // init load dictionary
  print_memory_usage("Memory before dictionary load");

  dict* d = dict_new(NSLOTS_MY_NODE);
  double timer = wtime();
  FILE* fp = fopen("data/digests/0", "r");
  load_file_to_dict(d, fp);
  timer = wtime() - timer;

  double nMB = get_file_size(fp) / ((double) 1000000);
  double nelm_sec = d->nelements_asked_to_be_inserted /  timer;

  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	 "dict read in %0.2fsec\n"
	 "It has %lu elms, file has %lu elms\n"
	 "d->nslots = %lu, d->nelements=%lu, filling rate=%f \n"
	 "i.e. it reads %0.4f MB/sec, and %0.4f elm /sec\n"
	 "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
	 timer,
	 d->nelements, d->nelements_asked_to_be_inserted,
	 d->nslots, d->nelements,
	 ((float) d->nelements)/d->nslots,
	 nMB/timer,
	 nelm_sec);


  print_memory_usage("After loadin dict");
  fclose(fp);
  

  // query large buffer
  u64 nqueries = PROCESS_QUOTA*10;
  int factor = 10;
  
  u64 nfound = 0;
  u8* queries = (u8*) malloc(sizeof(u8)*nqueries*N + 8);
  double elapsed = wtime();
  double start = 0;

  
  for (int j=0; j<10; j++) {
    getrandom(queries, nqueries*N, 0);

    start = wtime();
    for (size_t i=0; i<factor; i++)
     nfound += dict_has_elm(d, &queries[i*N]);

    elapsed += wtime() - start;

  }
  elapsed = wtime() - elapsed;
  double query_size_MB = N*nqueries / ((double) 1000000) ;
  
  printf("ignore me: nfound = %llu\n", nfound);
  printf("Querying %llu, took %0.2f sec i.e. %0.2f elm/sec = 2^%0.4f elm/sec \n"
	 "i.e %0.4f MB/sec\n",
	 nqueries*factor,
	 elapsed, nqueries*factor / elapsed, log2(nqueries*factor / elapsed),
	 query_size_MB / elapsed);
  
  
  
  

  return 0;

       
}

// example output:
/* Estimated memory per dictionary=2147483712 Bytes ≈ 2^31.00 Bytes */
/* thd2 sha_avx512_16way  elapsed 0.22sec i.e. 7187296.34 hashes/sec = 2^22.777 hashes, 459.9870 MB */
/* thd0 sha_avx512_16way  elapsed 0.22sec i.e. 7131176.55 hashes/sec = 2^22.766 hashes, 456.3953 MB */
/* thd1 sha_avx512_16way  elapsed 0.23sec i.e. 7013087.35 hashes/sec = 2^22.742 hashes, 448.8376 MB */
/* thd3 sha_avx512_16way  elapsed 0.22sec i.e. 7140152.36 hashes/sec = 2^22.768 hashes, 456.9698 MB */
/* thd5 sha_avx512_16way  elapsed 0.23sec i.e. 7009153.94 hashes/sec = 2^22.741 hashes, 448.5859 MB */
/* thd6 sha_avx512_16way  elapsed 0.22sec i.e. 7168525.04 hashes/sec = 2^22.773 hashes, 458.7856 MB */
/* thd7 sha_avx512_16way  elapsed 0.22sec i.e. 7204711.75 hashes/sec = 2^22.781 hashes, 461.1016 MB */
/* thd4 sha_avx512_16way  elapsed 0.23sec i.e. 7059654.51 hashes/sec = 2^22.751 hashes, 451.8179 MB */
/* sizeof(dict)=64 bytes */
/* Memory before dictionary loadMemory usage: ram: 801204 kb vm: 1334960 kb */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* dict read in 42.34sec */
/* It has 514064653 elms, file has 536871443 elms */
/* d->nslots = 536870912, d->nelements=514064653, filling rate=0.957520  */
/* i.e. it reads 139.4832 MB/sec, and 12680294.7343 elm /sec */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* After loadin dictMemory usage: ram: 2898508 kb vm: 5529268 kb */
/* ignore me: nfound = 0 */
/* Querying 100000000, took 2.07 sec i.e. 48196535.69 elm/sec = 2^25.5224 elm/sec  */
/* i.e 53.0162 MB/sec */

/* 256 Probes max */
/* Estimated memory per dictionary=2147483712 Bytes ≈ 2^31.00 Bytes */
/* sizeof(dict)=64 bytes */
/* Memory before dictionary loadMemory usage: ram: 1760 kb vm: 7328 kb */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* dict read in 46.24sec */
/* It has 525888347 elms, file has 536871443 elms */
/* d->nslots = 536870912, d->nelements=525888347, filling rate=0.979543  */
/* i.e. it reads 127.7057 MB/sec, and 11609607.3433 elm /sec */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* After loadin dictMemory usage: ram: 2099460 kb vm: 4201636 kb */
/* ignore me: nfound = 0 */
/* Querying 100000000, took 2.07 sec i.e. 48229017.18 elm/sec = 2^25.5234 elm/sec  */
/* i.e 53.0519 MB/sec */


/* 128 probes max */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* dict read in 43.89sec */
/* It has 521135438 elms, file has 536871443 elms */
/* d->nslots = 536870912, d->nelements=521135438, filling rate=0.970690  */
/* i.e. it reads 134.5424 MB/sec, and 12231126.3016 elm /sec */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* After loadin dictMemory usage: ram: 2099480 kb vm: 4201636 kb */
/* ignore me: nfound = 0 */
/* Querying 100000000, took 2.09 sec i.e. 47896585.59 elm/sec = 2^25.5134 elm/sec  */
/* i.e 52.6862 MB/sec */


/* 64 probes max */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* dict read in 42.48sec */
/* It has 514064653 elms, file has 536871443 elms */
/* d->nslots = 536870912, d->nelements=514064653, filling rate=0.957520  */
/* i.e. it reads 139.0044 MB/sec, and 12636766.7437 elm /sec */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* After loadin dictMemory usage: ram: 2099508 kb vm: 4201636 kb */
/* ignore me: nfound = 0 */
/* Querying 100000000, took 2.08 sec i.e. 48028335.81 elm/sec = 2^25.5174 elm/sec  */
/* i.e 52.8312 MB/sec */

/* 32 Probes max */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* dict read in 37.54sec */
/* It has 502976864 elms, file has 536871443 elms */
/* d->nslots = 536870912, d->nelements=502976864, filling rate=0.936867  */
/* i.e. it reads 157.3283 MB/sec, and 14302573.1757 elm /sec */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* After loadin dictMemory usage: ram: 2099464 kb vm: 4201636 kb */
/* ignore me: nfound = 0 */
/* Querying 100000000, took 2.07 sec i.e. 48363430.49 elm/sec = 2^25.5274 elm/sec  */
/* i.e 53.1998 MB/sec */
