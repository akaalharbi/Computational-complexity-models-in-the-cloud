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
#include "util_files.h"
#include <sys/random.h> // getrandom(void *buffer, size_t length, 1)
#include <mpi.h>
#include "common.h"
#include "sender.h"
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include "c_sha256_avx.h"
#include "arch_avx512_type1.h"



#define N_ACTIVE_LANES 16 /* this should be somewhere else */

inline static void increment_message(u8 Mavx[16][HASH_INPUT_SIZE])
{
  for (int i=0; i<(AVX_SIZE/WORD_SIZE_BITS); ++i) 
    ((CTR_TYPE*)Mavx[i])[0] += (i+1); /* increase counter part in M by 1 */
}

/* @todo move these function to util_arrays */
void print_m512i_u32(__m512i a, char* text){
  
  uint32_t A[16] = {0};
  _mm512_storeu_si512 ((__m256i*)A, a);
  printf("%s = ", text);
  for (int i = 0; i<16; ++i) {
    printf("%08x, ", A[i]);
  }
  puts("");
}

void print_2d_u32(u32** restrict a, int n, int m){
  for (int i=0; i<n; ++i) {
    for (int j=0; j<m; ++j) {
      printf("%x, ", a[i][j]);
    }
    puts("");
  }
}


void print_u32(u32* a, size_t l){
  for (size_t i = 0; i<l; ++i) 
    printf("%x, ", a[i]);
  puts("");
}



void extract_dist_points_dynamic(WORD_TYPE tr_states[restrict 16 * NWORDS_STATE],
			 int n_active_lanes, /* in */
			 u8 Mavx[restrict 16][HASH_INPUT_SIZE],
			 u8 digests[restrict 16 * N], /* out */
			 CTR_TYPE msg_ctrs_out[restrict 16], /* out */
			 int* n_dist_points) /* out */
{

  // ==========================================================================+
  // Summary: Find the distinguished digests in tr_state, and save them to     |
  //          states. Also, return how many disitingusihed points have been    |
  //          found. A point is disitingusihed if it has x bits on the right   |
  //          most that are zero. x = `DIFFICULTY` defined in config.h         |
  // INPUT:                                                                    |
  // - tr_states: transposed states of hash (specific  to sha256)              |
  // - max_nlanes: How many hashes out of 16 we consider. e.g. only first 7    | 
  // --------------------------------------------------------------------------+


  /* Init an AVX512/AVX2 vector */
  const REG_TYPE zero = SIMD_SETZERO_SI();
  /* use this mask to check if a digest is a distinguished point or not! */
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(DIST_PT_MASK);
  static REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  static REG_TYPE cmp_vect;
  static u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */
  int npotenital_dist_pt = 0;
  // test for distinguishedn point //

  /* load the last words of digests, we assume digest is aligned  */
  /* load last significant row in tr_state i.e. last words of each digest */
  digests_last_word = SIMD_LOAD_SI(&tr_states[(N_NWORDS_CEIL - 1)*16]);

  /* A distinguished point will have cmp_vect ith entry =  0  */
  cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);
  /* cmp_mask will have the ith bit = 1 if the ith element in cmp_vect is -1 */
  // Please the if is in one direction. 
  cmp_mask = SIMD_CMP_EPI32(cmp_vect, zero);
  /* printf("cmp_mask = %u\n", cmp_mask); */


  if (cmp_mask) { /* found at least a distinguished point? */
    npotenital_dist_pt = __builtin_popcount(cmp_mask);
    int lane = 0;
    int trailing_zeros = 0;

    for (int i=0;
	 (i < npotenital_dist_pt) && (lane < n_active_lanes);
	 ++i)
      {
      /* Basically get the index of the set bit in cmp_mask */
      /* Daniel Lemire has other methods that are found in his blog */
      trailing_zeros = __builtin_ctz(cmp_mask); 
      lane += trailing_zeros;
      cmp_mask = (cmp_mask >> trailing_zeros) ^ 1;

      if (lane < n_active_lanes){
	/* do usefule work only when we are acting on active lane  */
	*n_dist_points = i;
        /* update counter the ith counter */
	msg_ctrs_out[i] = ((CTR_TYPE*)Mavx[lane])[0];
	/* get the digest to digests vector */
	copy_transposed_digest(&digests[i*N], tr_states, lane);
      }
      
    }
  } /* end if (cmp_mask) */
}


void extract_dist_points(WORD_TYPE tr_states[restrict 16 * NWORDS_STATE],
			 u8 Mavx[restrict 16][HASH_INPUT_SIZE],
			 u8 digests[restrict 16 * N], /* out */
			 CTR_TYPE msg_ctrs_out[restrict 16], /* out */
			 int* n_dist_points) /* out */
{

  // ==========================================================================+
  // Summary: Find the distinguished digests in tr_state, and save them to     |
  //          states. Also, return how many disitingusihed points have been    |
  //          found. A point is disitingusihed if it has x bits on the right   |
  //          most that are zero. x = `DIFFICULTY` defined in config.h         |
  // INPUT:                                                                    |
  // - tr_states: transposed states of hash (specific  to sha256)              |
  // - max_nlanes: How many hashes out of 16 we consider. e.g. only first 7    | 
  // --------------------------------------------------------------------------+


  /* Init an AVX512/AVX2 vector */
  const REG_TYPE zero = SIMD_SETZERO_SI();
  /* use this mask to check if a digest is a distinguished point or not! */
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(DIST_PT_MASK);
  REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  REG_TYPE cmp_vect;
  u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */

  // test for distinguishedn point //

  /* load the last words of digests, we assume digest is aligned  */
  /* load last significant row in tr_state i.e. last words of each digest */
  digests_last_word = SIMD_LOAD_SI(&tr_states[(N_NWORDS_CEIL - 1)*16]);

  /* A distinguished point will have cmp_vect ith entry =  0  */
  cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);
  /* cmp_mask will have the ith bit = 1 if the ith element in cmp_vect is -1 */
  // Please the if is in one direction. 
  cmp_mask = SIMD_CMP_EPI32(cmp_vect, zero);
  /* printf("cmp_mask = %u\n", cmp_mask); */


  if (cmp_mask) { /* found at least a distinguished point? */
    *n_dist_points = __builtin_popcount(cmp_mask);
    int lane = 0;
    int trailing_zeros = 0;

    for (int i=0;
	 (i < *n_dist_points);
	 ++i)
      {
      /* Basically get the index of the set bit in cmp_mask */
      /* Daniel Lemire has other methods that are found in his blog */
      trailing_zeros = __builtin_ctz(cmp_mask); 
      lane += trailing_zeros;
      cmp_mask = (cmp_mask >> trailing_zeros) ^ 1;


      /* do usefule work only when we are acting on active lane  */

      /* update counter the ith counter */
      msg_ctrs_out[i] = ((CTR_TYPE*)Mavx[lane])[0];
      /* get the digest to digests vector */
      copy_transposed_digest(&digests[i*N], tr_states, lane);

      
    }
  } /* end if (cmp_mask) */
}





void find_hash_distinguished(u8 Mavx[restrict 16][HASH_INPUT_SIZE],
			     u8  digests[restrict 16 * N]/* out*/,
			     CTR_TYPE msg_ctrs[restrict 16], /* out */
			     CTR_TYPE* restrict  ctr /* in , out */,
			     int*  n_dist_points /* out */)

{

  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found. Mutates M  |
  //          and resests Mstate in the process. Suitable for phase ii use only|
  // --------------------------------------------------------------------------+
  // We start with message M, computed hash(M) if it is not a distinguished    |
  // point then we change M slightly. Increment the first 64 bits of M by 1    |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - Mavx[lane][HASH_INPUT_SIZE]: place holder for messages to be hashed     |
  //      using avx. they are copies of M, however counter part  increases by  |
  //      one for each lane.                                                   |
  // - digests[16*N]: Record the resulted digest of the distinguished points   |
  //  Record distinguished points in order. It contains at least 1 dist point  |
  // - ctr : update the last counter we found
  // --------------------------------------------------------------------------+
  // WARNING: this function can't deal with more than 32zeros as dist_test     |
  // NOTE: |
  // --------------------------------------------------------------------------+
  // TODO: use hash_multiple instead                                           |
  // --------------------------------------------------------------------------+

  /* no need to construct init state with each call of the function */ 
  static u32* tr_states; /* store the resulted digests here  */
  
  /* copy the counter part of M */

  *n_dist_points = 0;
  
  while (1) { /* loop till a dist pt found */

    /* */
    increment_message(Mavx);
    *ctr += (AVX_SIZE/WORD_SIZE_BITS); /* update counter */


    /* Hash multiple messages at once */
    #ifdef  __AVX512F__
    /* HASH 16 MESSAGES AT ONCE */
    tr_states = sha256_multiple_x16(Mavx);  
    #endif

    #ifndef  __AVX512F__
    #ifdef    __AVX2__
    /* HASH 16 MESSAGES AT ONCE */
    tr_states = sha256_multiple_oct(Mavx);
    #endif
    #endif

    /* Check if tr_states has  distinguished point(s) */
    extract_dist_points(tr_states, Mavx, digests, msg_ctrs, n_dist_points);
    if (*n_dist_points) {
      return; /* found a distinguished point */
    }
  } /* end while(1) */
}





static void regenerate_long_message_digests(u8 Mavx[restrict 16][HASH_INPUT_SIZE],
					    u32 tr_states[restrict 16*8],
					    u32 un_tr_states[restrict 16*8],
					    u8  digests[restrict 16 * N],
					    CTR_TYPE msg_ctrs[restrict 16],
					    u8* work_buf,
					    u8* snd_buf,
					    size_t servers_counters[restrict NSERVERS],
					    int myrank,
					    int nsenders,
					    MPI_Comm inter_comm)
{
  // ==========================================================================+
  // Summary: Regenrate hashes using the states file, read specific part of the|
  //          file that only this process will read. Send digests to severs.   |
  // --------------------------------------------------------------------------+
  /* M = 64bit ctr || 64bit nonce || random value */
  /* 16 messages for avx, each message differs om the other in the counter part */
  /* except counter, they are all the same */
  /* u8* bsnd_buf = (u8*) malloc((((sizeof(u8)*N + sizeof(CTR_TYPE))* PROCESS_QUOTA) + MPI_BSEND_OVERHEAD) */
  /* 			      * NSERVERS); */

  // -------------------------------- PART 1 ----------------------------------+
  // VARIABLES DEFINITIONS
  
  MPI_Status status;
  MPI_Request request;
  int first_time = 1; /* 1 if we have not sent anything yet, 0 otherwise */

  FILE* fp = fopen("data/states", "r");
  int server_id, n_dist_points;
  double elapsed = 0;
  
  /* u8 Mavx[16][HASH_INPUT_SIZE] = {0}; */
  /* u32 tr_states[16*8] = {0}; /\* same as current_states but transposed *\/ */

  // -------------------------------- PART 2 ----------------------------------+
  // Where do I read from the states file? How many states should I work on?
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = (myrank * nstates)/nsenders;
  size_t end = ((myrank + 1) * nstates)/nsenders;
  /* how many lanes are we using for hasing? sometimes we don't have enough*/
  int n_active_lanes = 16;  /* hashes to use all lanes */
  
  size_t global_idx; /* where are we in the states file */
  size_t local_idx; /* where are we in the buffer copied from states file  */
  size_t idx; /* */
  int inited = 0; /* 0 if we need to clear the avx register */
  

  if (myrank == (nsenders-1))
    end = nstates; /* get the rest of states */


  /* get all states that i should work on, pad with 0s to get 0 mod 16 */
  WORD_TYPE* states = (WORD_TYPE*) malloc( ((end - begin) + ((begin - end)%16))
					   * sizeof(WORD_TYPE)
					   * NWORDS_STATE);

  printf("pid=%d sender%d before states memset\n", getpid(), myrank);
  /* is it important to initialize it with zeros? */
  memset(states,
	 0,
	 ((end - begin) + ((begin - end)%16))
	 * sizeof(WORD_TYPE)
	 * NWORDS_STATE);
  
  /* only load states that i am going to work on */
  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, HASH_STATE_SIZE, (end - begin), fp);
  fclose(fp);

  printf("sender%d, begin=%lu, end=%lu, nstates=%lu\n",
	 myrank, begin, end, nstates);
  
  /* Attached buffered memory to MPI process  */

  // -------------------------------- PART 3 ----------------------------------+
  // Regenrate the long message: quota = (end - begin)
  // Part a:  quota = 16x + r, treat 16x states now, then work on r states later
  
  /* Hash the long message again, 16 at a time. The remaining will  */
  for (global_idx = begin; global_idx < end; global_idx += 16){
    /* local_idx = 0 -> (end-global)/16 */
    local_idx = global_idx - begin ;
    n_active_lanes = MIN((end - global_idx), 16);

    
    inited = 0; /* tell sha256_x16 to copy the state */



    /* Read fresh 16 states, and put them in transposed way  */
    /* if the n_active_lanes < 16, states has extra padded zeros */
    transpose_state(tr_states, &states[local_idx*NWORDS_STATE]);
    
    /* set message counters for all lanes (although not all will be used  ) */
    for (int lane = 0; lane<16; ++lane)
      ((u64*) Mavx[lane])[0] = INTERVAL * (global_idx + lane);

    elapsed = wtime(); /* how long does it take to hash an interval */
    
    for (u64 hash_n=0; hash_n <((u64) INTERVAL); ++hash_n){
      /* hash 16 messages and copy it to tr_states  */
      // todo fix me please 
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states, inited),
	     16*HASH_STATE_SIZE);

      inited = 1; /* sha256_x16 has alread a copy of the state */
      
      /* increment the message counters after hashing */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;




      /* we can avoid branching here by using extract_dist_points_dynamic */
      /* however in most cases we are not going to use it unless we are in */
      /* the boundary. */
      if (n_active_lanes == 16){
	extract_dist_points(tr_states, /* transposed states */
			    Mavx, /* messages used */
			    digests, /* save the distinguished hashes here */
			    msg_ctrs, /* messages are the same except counter */
			    &n_dist_points); /* how many dist points found? */

      }else {
	extract_dist_points_dynamic(tr_states, /* transposed states */
			    n_active_lanes,
			    Mavx, /* messages used */
			    digests, /* save the distinguished hashes here */
			    msg_ctrs, /* messages are the same except counter */
			    &n_dist_points); /* how many dist points found? */
      }


	
      /* put the distinguished points in specific serverss buffer */
      for (int i = 0; i<n_dist_points; ++i){
	/* skip this iteration if it is not a distinguish point */
	if (is_dist_state( (u8*) &un_tr_states[i*NWORDS_STATE]) )
	  continue;
	
	server_id = to_which_server(&digests[i*N]);

	/* where to put digest in server_id buffer? */
	idx = servers_counters[server_id];

        /* copy the digest to the server buffer */
        memcpy(&work_buf[server_id*PROCESS_QUOTA*N
			+ idx*N],
	       &digests[N*i],
	       N);

        ++servers_counters[server_id];

	/* if the server buffer is full send immediately */
	if (servers_counters[server_id] == PROCESS_QUOTA){
	  /* printf("before sender%d, global_idx=%lu, to %d \n", */
	  /* 	 myrank, */
	  /* 	 global_idx, */
	  /* 	 server_id); */
	  
	  /* only call wait when its not the first message */
	  if (!first_time) 
	    MPI_Wait(&request, &status);

	  /* copy the message to be sent to snd_buf */
	  memcpy(snd_buf,
		 &work_buf[server_id*PROCESS_QUOTA*N], /* idx of server buf */
		 N*PROCESS_QUOTA);



	  MPI_Isend(snd_buf, 
		    N*PROCESS_QUOTA, /* How many characteres will be sent. */
		    MPI_UNSIGNED_CHAR, 
		    server_id, /* receiver */
		    TAG_DICT_SND, /* 0 */
		    inter_comm,
		    &request);

	  first_time = 0; /* we have sent a message */
	  
	  servers_counters[server_id] = 0;
	  
	}/* end if */

	
      } /* end for n_dist_points */
    } /* end for hash_n */

    elapsed = wtime() - elapsed;
    printf("sender%d, global_idx=%lu, 2^%f hashes/sec, done %f%%, %0.4fsec\n",
	   myrank,
	   global_idx,
	   log2(INTERVAL/elapsed)+log2(16),
	   100 * ((double_t) 1 -  (end - global_idx)/((double_t) end - begin)),
	    elapsed);
    
    elapsed = wtime();
    
  } /* end for global_idx */

  /* should we wait for the last message to be sent?! */
  if (!first_time) 
    MPI_Wait(&request, &status);


  /* Tell every receiver that you are done!, and send the remaining digests */
  for (int server=0; server<NSERVERS; ++server){
    /* we may have some digests that we have not sent since their number is */
    /* less than the PROCESS_QUOTA. We send them */

    // Before:
    /* server i buffer: {dgst1, .., dgstk, random garbage}*/
    /* lft term: server i buffer begin, right term: nbytes to reach garbage */
    memset(&work_buf[N*PROCESS_QUOTA*server + N*servers_counters[server]],
	   0,
	   N*(PROCESS_QUOTA-servers_counters[server]));
    // After:
    /* server i buffer: {dgst1, .., dgstk, 0, ..., 0}, len = PROCESS_QUOTA */
    /* By default the receiver will ignore zero digests */

    
    MPI_Send(&work_buf[server*PROCESS_QUOTA*N],
	     PROCESS_QUOTA*N,
	     MPI_UNSIGNED_CHAR,
	     server,
	     TAG_DONE_HASHING,
	     inter_comm);

    printf("sender%d dit au revoir\n", myrank);

  }

  /* أترك المكان كما كان أو أفضل ما كان  */
  memset(work_buf, 0, N*PROCESS_QUOTA*NSERVERS); 
  
}



static void generate_random_digests(u8 Mavx[16][HASH_INPUT_SIZE],/* random msg */
				    u8  digests[restrict 16 * N],
				    CTR_TYPE msg_ctrs[restrict 16],
				    u8* restrict work_buf,
				    u8* bsnd_buf,
				    size_t servers_counters[restrict NSERVERS],
				    MPI_Comm inter_comm)
{
  u32* states_avx;
  int n_dist_points = 0; /* At the beginning we no dist points */
  size_t server_id=0, idx=0;

  MPI_Status status;
  MPI_Request request;

  int first_send = 1;

  


  while (1) {
    /* 1- hash, 2- increment message counters, 3- extract disit point if any */
    states_avx = sha256_multiple_x16(Mavx);
    increment_message(Mavx);
    extract_dist_points(states_avx, Mavx, digests, msg_ctrs, &n_dist_points);

    /* put the distinguished points in specific serverss buffer */
    for (int i = 0; i<n_dist_points; ++i){ /* n_dist_points might be 0 */
      server_id = to_which_server(&digests[i*N]);
      /* where to put digest in server_id buffer? (this is a local view) */
      idx = servers_counters[server_id];


      /* 1 element = ctr||dgst  */
      /* #elments in  server buffer = PROCESS_QUOTA */
      /* There are NSERVERS in total */
      /* copy the ctr to the server buffer */
      memcpy(&work_buf[( N + sizeof(CTR_TYPE)  )/* one element size */
		       *( PROCESS_QUOTA*server_id + idx)],/* #elments skipped */
	     &msg_ctrs[i], /* counter of the message i  */
	     sizeof(CTR_TYPE));

      /* copy the digest to the server buffer */
      memcpy(&work_buf[( N + sizeof(CTR_TYPE)  ) /* one element size */
		       *( PROCESS_QUOTA*server_id + idx) /* #elments skipped */
		       + sizeof(CTR_TYPE)], /* write after the counter part */
	     &digests[N*i],
	     N);

      ++servers_counters[server_id];
      /* if the server buffer is full send immediately */
      if (servers_counters[server_id] == PROCESS_QUOTA){
	if (!first_send)
	  MPI_Wait(&request, &status);
	
	MPI_Isend(&work_buf[server_id*PROCESS_QUOTA*(N+sizeof(CTR_TYPE))],
		  (N+sizeof(CTR_TYPE))*PROCESS_QUOTA, /* #chars to be sent */
		  MPI_UNSIGNED_CHAR, 
		  server_id, /* receiver */
		  TAG_DICT_SND, /* 0 */
		  inter_comm,
		  &request);
	first_send = 0; /* call mpi wait next time */
	servers_counters[server_id] = 0;
      }
    }
  }
  
}



void sender( MPI_Comm local_comm, MPI_Comm inter_comm)
{

  // ------------------------------------------------------------------------+
  // I am a sending processor. I only generate hashes and send them.         |
  // Process Numbers: [NSERVERS + 1,  nproc]                                 |
  //-------------------------------------------------------------------------+

  // ----------------------------- PART 0 --------------------------------- //
  int myrank, nsenders;
  MPI_Comm_rank(local_comm, &myrank);
  MPI_Comm_size(local_comm, &nsenders);

  /* random message base */
  u8 M[HASH_INPUT_SIZE];
  /* non-transposed messages */
  u8 Mavx[16][HASH_INPUT_SIZE] = {0};
  /* transoposed states  */
  u32 tr_states[16*8] __attribute__ ((aligned(64))) = {0};
  u32 un_tr_states[16*8] __attribute__ ((aligned(64))) = {0};  /* non-transoposed */
  u8 digests[16 * N] ; /* save the distinguished digests here (manually) */
  CTR_TYPE msg_ctrs[16]; /* save thde counters of the dist digests (manually) */

  

  /* Sending related buffers allocation  */

  /* |dgst| + |ctr| */
  size_t const one_pair_size = sizeof(u8)*N + sizeof(CTR_TYPE); 
  size_t const  nbytes_per_server = one_pair_size * PROCESS_QUOTA;

  /* pair = (ctr||dgst) */
  /* { (server0 paris..) | (server1 pairs...) | ... | (serverK pairs...) } */
  u8* work_buf = (u8*) malloc(nbytes_per_server * NSERVERS );

  /* attach to MPI_BSend, not worked directly on! */
  u8* bsnd_buf = (u8*) malloc((nbytes_per_server + MPI_BSEND_OVERHEAD)
			      * NSERVERS);
  /* ith_entry : how many non-sent element stored in server i buffer */
  size_t* server_counters = malloc(sizeof(size_t)*NSERVERS);



  if (bsnd_buf == NULL)
    puts("bsnd_buf is NULL");

  if (work_buf == NULL)
    puts("work_buf is NULL");

  printf("pid=%d sender%d before memset\n", getpid(), myrank);
  /* memset(bsnd_buf, 0, (nbytes_per_server + MPI_BSEND_OVERHEAD)* NSERVERS); */
  memset(work_buf, 0, nbytes_per_server*NSERVERS);
  memset(server_counters, 0, sizeof(size_t)*NSERVERS);
  

  printf("sender%d: N=%d, one pair size = %lu\n",
	 myrank, N, one_pair_size);

  
  // ----------------------------- PART 1 --------------------------------- //
  // Regenrate the long message in parallel!                                //
  regenerate_long_message_digests(Mavx,
				  tr_states,
				  un_tr_states,
				  digests,
				  msg_ctrs,
				  work_buf,
				  bsnd_buf,
				  server_counters,
				  myrank,
				  nsenders,
				  inter_comm);


  // ----------------------------- PART 2.a ----------------------------------- //
  /* Get a random message only once */
  CTR_TYPE*ctr_pt = (CTR_TYPE*) M; /* counter pointer */
  getrandom(M, HASH_INPUT_SIZE, 1);
  ctr_pt[0] = 0; /* zeroing the first 64bits of M */

  // copy the the random message to all avx messages
  for (int i = 0; i<16; ++i)
    memcpy(Mavx[i], M, HASH_INPUT_SIZE);

  // print the template. this is not necessary.
  char txt[50];
  snprintf(txt, sizeof(txt), "sender #%d template=", myrank);
  print_byte_txt(txt, M,HASH_INPUT_SIZE);
  puts("\n");

  
  // ----------------------------- PART 2.b ---------------------------------- //
  // 1-  Sen the initial input to all receiving servers
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(M, HASH_INPUT_SIZE, MPI_UNSIGNED_CHAR, NULL, 0, MPI_UNSIGNED_CHAR, inter_comm);

  printf("sender%d done sending template\n", myrank);

  // ------------------------------ PART 3 ----------------------------------- +
  // 1- generate disitingusihed hashes                                         |
  // 2- decide which server should probe the disitingusihed hash               |
  // 3- if server i buffer has enough hashes, send them immediately.           |
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}                               |
  // dgst := N bytes of the digest                                             |
  // ctr  := (this is a shortcut that allows us not to send the whole message) |
  //---------------------------------------------------------------------------+
  generate_random_digests(Mavx,
			  digests,
			  msg_ctrs,
			  work_buf,
			  bsnd_buf,
			  server_counters,
			  inter_comm);



  free(work_buf);
  free(bsnd_buf);
  free(server_counters);
  free(ctr_pt);
  

  return; // au revoir.

} // MPI_Finalize


