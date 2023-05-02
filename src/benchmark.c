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







int main(int argc, char *argv[])
{





  /* Benchmark dictionary query  */
  // init load dictionary

  size_t factor = 1000;
  size_t dict_size = PROCESS_QUOTA*factor;
  dict* d = dict_new(dict_size);
  dict* d_simd = dict_new(dict_size);
  u64 ndigests = PROCESS_QUOTA;
  u64 size_digests = ndigests*N;
  u8* message = malloc(size_digests);

  printf("d->nslots=%lu, d_simd->nslost=%lu\n", d->nslots, d_simd->nslots);


  double timer = 0;
  double elapsed_dict = 0;


  for (int i=0; i<factor; i++) {
    getrandom(message, size_digests, 0);
    timer = wtime();
    for (size_t i=0; i<PROCESS_QUOTA; ++i) {
      dict_add_element_to(d, (u8*) &message[N*i]);
    }
    elapsed_dict += wtime() - timer;

  }
  getrandom(message, PROCESS_QUOTA*N, 0);


  printf("Dictionary addtion took %0.2f i.e. 2^%0.4felm/sec\n",
	 elapsed_dict,
	 log2(dict_size/elapsed_dict));



  
  /* FILE* fp = fopen("data/digests/0", "r"); */
  /* load_file_to_dict(d, fp); */
  /* timer = wtime() - timer; */

  /* double nMB = get_file_size(fp) / ((double) 1000000); */
  /* double nelm_sec = d->nelements_asked_to_be_inserted /  timer; */

  /* printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" */
  /* 	 "dict read in %0.2fsec\n" */
  /* 	 "It has %lu elms, file has %lu elms\n" */
  /* 	 "d->nslots = %lu, d->nelements=%lu, filling rate=%f \n" */
  /* 	 "i.e. it reads %0.4f MB/sec, and %0.4f elm /sec\n" */
  /* 	 "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", */
  /* 	 timer, */
  /* 	 d->nelements, d->nelements_asked_to_be_inserted, */
  /* 	 d->nslots, d->nelements, */
  /* 	 ((float) d->nelements)/d->nslots, */
  /* 	 nMB/timer, */
  /* 	 nelm_sec); */


  /* print_memory_usage("After loadin dict"); */
  /* fclose(fp); */
  

  // query large buffer
  u64 nqueries = PROCESS_QUOTA*10;
  factor = 10;
  
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
