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
#include "c_sha256_avx.h"
#include "arch_avx512_type1.h"

inline static void increment_message(u8 Mavx[16][HASH_INPUT_SIZE])
{
  for (int i=0; i<(AVX_SIZE/WORD_SIZE_BITS); ++i) 
    ((CTR_TYPE*)Mavx[i])[0] += (i+1); /* increase counter part in M by 1 */
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
  // --------------------------------------------------------------------------+


  /* Init an AVX512/AVX2 vector */
  const REG_TYPE zero = SIMD_SETZERO_SI();
  /* use this mask to check if a digest is a distinguished point or not! */
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(MASK);
  static REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  static REG_TYPE cmp_vect;
  static u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */

  // test for distinguishedn point //

  /* load the last words of digests, we assume digest is aligned  */
  digests_last_word = SIMD_LOAD_SI(&tr_states[(N_NWORDS_CEIL - 1) * HASH_STATE_SIZE]);
  /* A distinguished point will have cmp_vect ith entry =  0  */
  cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);
  /* cmp_mask will have the ith bit = 1 if the ith element in cmp_vect is -1 */
  // Please the if is in one direction. 
  cmp_mask = SIMD_CMP_EPI32(cmp_vect, zero);

  if (cmp_mask) { /* found at least a distinguished point? */
    *n_dist_points = __builtin_popcount(cmp_mask);
    int lane = 0;
    int trailing_zeros = 0;

    for (int i=0; i < (*n_dist_points); ++i){
      /* Basically get the index of the set bit in cmp_mask */
      /* Daniel Lemire has other methods that are found in his blog */
      trailing_zeros = __builtin_ctz(cmp_mask); 
      lane += trailing_zeros;
      cmp_mask = (cmp_mask >> trailing_zeros) ^ 1;
       
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


static void regenerate_long_message_digests(u8 Mavx[16][HASH_INPUT_SIZE],
					    u8  digests[restrict 16 * N],
					    CTR_TYPE msg_ctrs[restrict 16],
					    u8* restrict work_buf,
					    u8* bsnd_buf,
					    size_t servers_counters[restrict NSERVERS],
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

  FILE* fp = fopen("data/states", "r");
  int myrank;
  size_t server_id, idx;
  int n_dist_points = 0;
  
  MPI_Comm_rank(inter_comm, &myrank);
  
  size_t nstates = get_file_size(fp) / HASH_STATE_SIZE;
  size_t begin = myrank * (nstates/nsenders);
  size_t end = (myrank + 1) * (nstates/nsenders);

  
  if (myrank == (nsenders-1))
    end = nstates;

  /* get all states that i should work on: */
  WORD_TYPE* states = (WORD_TYPE*) malloc((end - begin) * sizeof(WORD_TYPE));  
  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, WORD_SIZE, (end - begin), fp);

  /* Hash buffers init */
  SHA256_ARGS args;
  MPI_Buffer_attach(bsnd_buf,
		    PROCESS_QUOTA*N*NSERVERS
		    + NSERVERS*MPI_BSEND_OVERHEAD);
  
  /* Hash the long message again, 16 at a time */
  for (size_t step = 0; step<(end - begin)/16; ++step){
    /* Read fresh states, and put them in transposed way  */
    for (int lane = 0; lane < 16; lane++) {
      args.digest[lane + 0*16] = states[NWORDS_STATE*(step*16 + lane) + 0];
      args.digest[lane + 1*16] = states[NWORDS_STATE*(step*16 + lane) + 1];
      args.digest[lane + 2*16] = states[NWORDS_STATE*(step*16 + lane) + 2];
      args.digest[lane + 3*16] = states[NWORDS_STATE*(step*16 + lane) + 3];
      args.digest[lane + 4*16] = states[NWORDS_STATE*(step*16 + lane) + 4];
      args.digest[lane + 5*16] = states[NWORDS_STATE*(step*16 + lane) + 5];
      args.digest[lane + 6*16] = states[NWORDS_STATE*(step*16 + lane) + 6];
      args.digest[lane + 7*16] = states[NWORDS_STATE*(step*16 + lane) + 7];

      /* adjust the counter part */
      ((u64*) Mavx[lane])[0] = (step*16 + lane)*INTERVAL;
      args.data_ptr[lane] = Mavx[lane]; /* point data */
      
    }


    for (size_t hash_n=0; hash_n < INTERVAL; ++hash_n){
      /* hash 16 messages  */
      call_sha256_x16_avx512_from_c(&args, 1);
      /* update message counters */

      extract_dist_points(args.digest, /* transposed state */
			  Mavx, /* messages used */
			  digests, /* save the distinguished hashes here */
			  msg_ctrs, /* messages are the same except counter */
			  &n_dist_points); /* how many dist points found? */

      /* put the distinguished points in specific serverss buffer */
      for (int i = 0; i<n_dist_points; ++i){
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
	  MPI_Bsend(&work_buf[server_id*PROCESS_QUOTA*N], /* idx of server buf */
		    N*PROCESS_QUOTA, /* How many characteres will be sent. */
		    MPI_UNSIGNED_CHAR, 
		    server_id, /* receiver */
		    TAG_DICT_SND, /* 0 */
		    inter_comm);

	  servers_counters[server_id] = 0;
	}
      }
      
      
      increment_message(Mavx); /* changes should carry to args.data_ptr */
      /* for (int lane = 0; lane<16; ++lane){ */
      /* 	((u64*) Mavx[lane])[0] += 1; */
      /* 	args.data_ptr[lane] = Mavx[lane]; /\* not sure this is necessary *\/ */
      /* } */
      
    }
  }

  int attached_size;
  MPI_Buffer_detach(&bsnd_buf, &attached_size);
  /* أترك المكان كما كان أو أفضل ما كان  */
  memset(work_buf, 0, N*PROCESS_QUOTA*NSERVERS); 

  /* Tell every receiver that you are done! */
  for (int server=0; server<NSERVERS; ++server) 
    MPI_Send(work_buf,
	     PROCESS_QUOTA*N,
	     MPI_UNSIGNED_CHAR,
	     server,
	     TAG_DONE_HASHING,
	     inter_comm);
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
  int n_dist_points = 0; /* At the beginning we no dist points
			  */
  size_t server_id=0, idx=0;
  
  MPI_Buffer_attach(bsnd_buf,
		    PROCESS_QUOTA*(N+sizeof(CTR_TYPE))*NSERVERS
		    + NSERVERS*MPI_BSEND_OVERHEAD);


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
	  
	MPI_Bsend(&work_buf[server_id*PROCESS_QUOTA*(N+sizeof(CTR_TYPE))],
		  (N+sizeof(CTR_TYPE))*PROCESS_QUOTA, /* #chars to be sent */
		  MPI_UNSIGNED_CHAR, 
		  server_id, /* receiver */
		  TAG_DICT_SND, /* 0 */
		  inter_comm);
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
  u8 Mavx[16][HASH_INPUT_SIZE] = {0}; /* non-transposed! */
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
  u8* bsnd_buf = (u8*) malloc(nbytes_per_server * NSERVERS);
  /* ith_entry : how many non-sent element stored in server i buffer */
  size_t* server_counters = malloc(sizeof(size_t)*NSERVERS);



  

  
  // ----------------------------- PART 1 --------------------------------- //
  // Regenrate the long message in parallel!                                //
  regenerate_long_message_digests(Mavx,
				  digests,
				  msg_ctrs,
				  work_buf,
				  bsnd_buf,
				  server_counters,
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
  MPI_Allgather(M, HASH_INPUT_SIZE, MPI_UNSIGNED_CHAR, NULL, 0, NULL, inter_comm);


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


