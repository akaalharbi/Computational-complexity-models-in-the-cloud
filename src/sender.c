#include "numbers_shorthands.h"
#include "hash.h"
#include "sha256_single.h"
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


static inline void show_and_save_benchmark
    (double total_elapsed,
     double elapsed_mpi_send,
     double elapsed_hash,
     double elapsed_extract_dist,
     size_t nmsgs_sent,
     size_t nbytes_per_server,
     size_t nhashes,
     size_t interval,
     int nsenders,
     int myrank,
     FILE* fp)
{


  printf("->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->\n"
	 "total=%fsec, mpi_wait=%fsec, hash=%fsec≈%fhash/sec≈%fMB/sec, find dist=%fsec\n"
	 "mpi_send=%f%%, hash=%f%%, find dist=%f%%\n" 
	 "send %f MB/sec, exp[all senders] = %f MB/sec, nsenders=%d, nservers=%d\n"
	 "DIFFICULTY=%d, INTERVAL=%d, nsends=%lu\n"
	 "->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->\n",
	 total_elapsed,
	 elapsed_mpi_send,
	 elapsed_hash,
	 (nhashes) / elapsed_hash,
	 (nhashes*N) / (elapsed_hash*1000000),
	 elapsed_extract_dist,
	 100*elapsed_mpi_send/total_elapsed,
	 100*elapsed_hash/total_elapsed,
	 100*elapsed_extract_dist/total_elapsed,
	 nmsgs_sent*((nbytes_per_server)/total_elapsed)/1000000,
	 nsenders*nmsgs_sent*((nbytes_per_server)/total_elapsed)/1000000,
	 nsenders,
	 NSERVERS,
	 DIFFICULTY,
	 (int) log2(interval),
	 nmsgs_sent);


  fprintf(fp,"->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->\n"
	 "total=%fsec, mpi_wait=%fsec, hash=%fsec≈%fhash/sec≈%fMB/sec, find dist=%fsec\n"
	 "mpi_send=%f%%, hash=%f%%, find dist=%f%%\n" 
	 "send %f MB/sec, exp[all senders] = %f MB/sec, nsenders=%d, nservers=%d\n"
	 "DIFFICULTY=%d, INTERVAL=%d, nsends=%lu\n"
	 "->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->\n",
	 total_elapsed,
	 elapsed_mpi_send,
	 elapsed_hash,
	 (nhashes) / elapsed_hash,
	 (nhashes*N) / (elapsed_hash*1000000),
	 elapsed_extract_dist,
	 100*elapsed_mpi_send/total_elapsed,
	 100*elapsed_hash/total_elapsed,
	 100*elapsed_extract_dist/total_elapsed,
	 nmsgs_sent*((nbytes_per_server)/total_elapsed)/1000000,
	 nsenders*nmsgs_sent*((nbytes_per_server)/total_elapsed)/1000000,
	 nsenders,
	 NSERVERS,
	 DIFFICULTY,
	 (int) log2(interval),
	 nmsgs_sent);

  
}    


void buffered_isend(u8 *snd_bufs, /* all send bufs */
		    MPI_Request requests[],
		    MPI_Status statuses[],
		    int nbufs, /* how many asynchrous send at the same time */    
		    u8* msg, /* msg to be sent */
		    int msg_length, /* length in bytes */
		    int to,
		    int TAG,
		    MPI_Comm comm)
{
  /* a simple implementation of personal view of how a buffered send work */
  // we have nbufs buffers that will be used to send asynchrously.
  // We will find the first one whose sending job has been completed
  // then we initiate a new send.
  
  int is_free = 0;
  int idx = -1;



  
  /* find a buffer with a completed sending job  */
  while (!is_free)
    MPI_Testany(nbufs, requests, &idx, &is_free, MPI_STATUS_IGNORE);
    
  idx = (idx >= 0) ? idx : 0; /* see https://github.com/open-mpi/ompi/issues/11725 */
  /* copy the message to be sent to the free send buffer */
  memcpy(&snd_bufs[idx*msg_length], msg, msg_length);

  MPI_Isend(&snd_bufs[idx*msg_length],
	    msg_length,
	    MPI_UNSIGNED_CHAR,
	    to,
	    TAG,
	    comm,
	    &requests[idx]);
  

  
    
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
  //          most position that are zero. x = `DIFFICULTY` defined in config.h|
  // INPUT:                                                                    |
  // - tr_states: transposed states of hash (specific  to sha256)              |
  // - n_active_lanes: How many hashes out of 16 we consider. e.g. only first 7|
  // - Mavx: Messages that used to generate this tr_states.                    |
  // - msg_ctrs_out: save the messages counters that when are hashed generate  |
  //                 distinguished points.                                     |
  // - n_dist_points: how many distinguished points have been found.           |
  // --------------------------------------------------------------------------+


  /* Init an AVX512 vector */
  const REG_TYPE zero = SIMD_SETZERO_SI();

  /* use this mask to check if a digest is a distinguished point or not! */
  /* digest & dist_pt_mask == 0 iff digest is distinguished point */
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(DIST_PT_MASK);
  REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  REG_TYPE cmp_vect;
  u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */
  /* since we restrict our treatment to only n_active_lane */
  int npotenital_dist_pt = 0;

  /* INIT */
  *n_dist_points = 0; /* old data should not interfer with new one */
  
  /* load last significant row in tr_state i.e. last words of each digest */
  digests_last_word = SIMD_LOADU_SI(&tr_states[(N_NWORDS_CEIL - 1)*16]);

  /* A distinguished point will have ith entry in  cmp_vect =  0  */
  cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);

  /* cmp_mask will have the ith bit = 1 if the ith element in cmp_vect is -1 */
  cmp_mask = SIMD_CMP_EPI32(cmp_vect, zero);

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
      /* remove all bits up to the found lane  */
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


void extract_dist_points_16(WORD_TYPE tr_states[restrict 16 * NWORDS_STATE],
			 u8 Mavx[restrict 16][HASH_INPUT_SIZE],
			 u8 digests[restrict 16 * N], /* out */
			 CTR_TYPE msg_ctrs_out[restrict 16], /* out */
			 int* n_dist_points) /* out */
{
  // ==========================================================================+
  // Summary: Find the distinguished digests in tr_state, and save them to     |
  //          states. Also, return how many disitingusihed points have been    |
  //          found. A point is disitingusihed if it has x bits on the right   |
  //          most position that are zero. x = `DIFFICULTY` defined in config.h|
  // INPUT:                                                                    |
  // - tr_states: transposed states of hash (specific  to sha256)              |
  // - n_active_lanes: How many hashes out of 16 we consider. e.g. only first 7|
  // - Mavx: Messages that used to generate this tr_states.                    |
  // - msg_ctrs_out: save the messages counters that when are hashed generate  |
  //                 distinguished points.                                     |
  // - n_dist_points: how many distinguished points have been found.           |
  // --------------------------------------------------------------------------+


  /* Init an AVX512/AVX2 vector */
  const REG_TYPE zero = SIMD_SETZERO_SI();

  /* use this mask to check if a digest is a distinguished point or not! */
  const REG_TYPE dist_mask_vect = SIMD_SET1_EPI32(DIST_PT_MASK);
  REG_TYPE digests_last_word ; /* _mm512_load_epi32 */
  REG_TYPE cmp_vect;
  u16 cmp_mask = 0; /* bit i is set iff the ith digest is disitingusihed */
  *n_dist_points = 0; /* old data should not interfer with the new one */
  
  // test for distinguishedn point //
  /* load the last words of digests, we assume digest is aligned  */
  /* load last significant row in tr_state i.e. last words of each digest */
  digests_last_word = SIMD_LOADU_SI(&tr_states[(N_NWORDS_CEIL - 1)*16]);

  /* A distinguished point will have cmp_vect ith entry =  0  */
  cmp_vect = SIMD_AND_EPI32(digests_last_word, dist_mask_vect);
  /* cmp_mask will have the ith bit = 1 if the ith element in cmp_vect is -1 */

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

      /* remove all bits up to the found lane  */
      cmp_mask = (cmp_mask >> trailing_zeros) ^ 1;

      /* update counter the ith counter */
      msg_ctrs_out[i] = ((CTR_TYPE*)Mavx[lane])[0];

      /* get the digest to digests vector */
      copy_transposed_digest(&digests[i*N], tr_states, lane);
      /* assert(is_dist_digest(&digests[i*N])); */
    }
  } /* end if (cmp_mask) */
}




static void regenerate_long_message_digests(u8 Mavx[restrict 16][HASH_INPUT_SIZE],
					    u32 tr_states[restrict 16*8],
					    u8  digests[restrict 16 * N],
					    CTR_TYPE msg_ctrs[restrict 16],
					    u8* work_buf,
					    u8* snd_bufs,
					    MPI_Request requests[],
					    MPI_Status statuses[],
					    int nbufs, /* #async sends */
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


  FILE* fp = fopen("data/states", "r");

  char file_name[FILE_NAME_MAX_LENGTH];
  snprintf(file_name, sizeof(file_name), "data/stats/sender_%d", myrank);
  FILE* fp_timing = fopen(file_name, "w");

  
  int server_id, n_dist_points;
  
  /* timing variables for profiliing sending, hashing, and extracting dist pt */
  double elapsed = 0; 
  double total_elapsed = 0;
  double timer = 0; /* general timer start */
  double elapsed_hash = 0;
  double elapsed_extract_dist = 0;
  double elapsed_mpi_send = 0;
  size_t nmsgs_sent = 0;
  
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

   
  /* only load states that i am going to work on */

  fseek(fp, begin*HASH_STATE_SIZE, SEEK_SET);
  fread(states, HASH_STATE_SIZE, (end - begin), fp);
  /* CLOSE THE FILE! */
  fclose(fp);

  printf("sender%d, begin=%lu, end=%lu, quotua=%lu, nstates=%lu\n",
	 myrank, begin, end, (end-begin), nstates);





  // -------------------------------- PART 3 ----------------------------------+
  // Regenrate the long message: quota = (end - begin)
  // Part a:  quota = 16x + r, treat 16x states now, then work on r states later
  total_elapsed = wtime();
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
      
      timer = wtime();
      memcpy(tr_states,
	     sha256_multiple_x16_tr(Mavx, tr_states, inited),
	     16*HASH_STATE_SIZE);
      elapsed_hash += (wtime() - timer);


      
      inited = 1; /* sha256_x16 has alread a copy of the state */
      
      /* increment the message counters after hashing */
      for (int lane = 0; lane<16; ++lane)
	((u64*) Mavx[lane])[0] += 1;


      /* we can avoid branching here by using extract_dist_points_dynamic */
      /* however in most cases we are not going to use it unless we are in */
      /* the boundary. */
      timer = wtime();
      if (n_active_lanes == 16){
	extract_dist_points_16(tr_states, /* transposed states */
			    Mavx, /* messages used */
			    digests, /* save the distinguished hashes here */
			    msg_ctrs, /* messages are the same except counter */
			    &n_dist_points); /* how many dist points found? */

      }else {/* this function slower than the one above! */
	extract_dist_points_dynamic(tr_states, /* transposed states */
			    n_active_lanes,
			    Mavx, /* messages used */
			    digests, /* save the distinguished hashes here */
			    msg_ctrs, /* messages are the same except counter */
			    &n_dist_points); /* how many dist points found? */
      }
      elapsed_extract_dist += wtime() - timer;

      
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

	/* if the server buffer is full, send immediately */
	if (servers_counters[server_id] == PROCESS_QUOTA){

	  timer = wtime();
	  /* send  &work_buf[server_id*PROCESS_QUOTA*N] -> server[server_id] */
	  buffered_isend(snd_bufs,
			 requests,
			 statuses,
			 nbufs,
			 &work_buf[server_id*PROCESS_QUOTA*N],
			 N*PROCESS_QUOTA,
			 server_id,
			 TAG_DICT_SND,
			 inter_comm);
	  
	  /* it should be called elapsed_mpi_send */
	  elapsed_mpi_send += wtime() - timer;

	  ++nmsgs_sent;
	  /* consider it as empty buffer without the overhead of erasing it */
	  servers_counters[server_id] = 0;
	}/* end if server_counters */
      } /* end for n_dist_points */
    } /* end for hash_n */

    elapsed = wtime() - elapsed;
    printf("sender%d, global_idx=%lu, 2^%f hashes/sec, done %0.4f%%, %0.4fsec, ETA %f sec\n",
	   myrank,
	   global_idx,
	   log2(INTERVAL/elapsed)+log2(16),
	   100 * ((double_t) 1 -  (end - global_idx)/((double_t) end - begin)),
	   elapsed,
	   elapsed * ( (end - global_idx)/16));
    elapsed = wtime();
  } /* end for global_idx */
  /* we have a hanging MPI_Isend, make sure it has been sent */
  
  total_elapsed = wtime() - total_elapsed;

  show_and_save_benchmark(total_elapsed,
			  elapsed_mpi_send,
			  elapsed_hash,
			  elapsed_extract_dist,
			  nmsgs_sent,
			  (N*PROCESS_QUOTA),
			  (end-begin)*INTERVAL,
			  INTERVAL,
			  nsenders,
			  myrank,
			  fp_timing);
  


  // -------------------------------- PART 4 ----------------------------------+
  // Send the digests that remained in the buffer and tell the receivers that
  // I finished hashing.

  /* check there are no more mpi_isend hanging */
  for (int k = 0; k<nbufs; ++k)
    MPI_Wait(&requests[k], &statuses[k]);


  /* Tell every receiver that you are done!, and send the remaining digests */
  for (int server=0; server<NSERVERS; ++server){
    /* we may have some digests that we have not sent since their number is */
    /* less than the PROCESS_QUOTA. We have to send them!.  */

    /* set all the empty digests places in work_buffer as zero since they will  */
    /* be ignored by the receiver by default. */
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
  }
  printf("sender%d dit au revoir\n", myrank);
  /* أترك المكان كما كان أو أفضل ما كان  */
  memset(work_buf, 0, N*PROCESS_QUOTA*NSERVERS); 
  memset(servers_counters, 0, sizeof(size_t)*NSERVERS);
  fclose(fp_timing);
}



static void generate_random_digests(u8 Mavx[16][HASH_INPUT_SIZE],/* random msg */
				    u8  digests[restrict 16 * N],
				    CTR_TYPE msg_ctrs[restrict 16],
				    u8* restrict work_buf,
				    u8* snd_bufs,
				    MPI_Request requests[],
				    MPI_Status statuses[],
				    int nbufs,
				    size_t servers_counters[restrict NSERVERS],
				    const size_t nbytes_per_server,
				    const size_t one_elm_size,
				    int myrank,
				    int nsenders,
				    MPI_Comm inter_comm)
{
  /// 1- increment message counters, 2- hash.  3- extract disit point if any,
  /// repeate. IMPORTANT: DO NOT CHANGE THE ORDER OF THIS SEQUENCE!

  /* timing variables for profiliing sending, hashing, and extracting dist pt */
  double total_elapsed = wtime(); /* save the current time */
  double timer = 0; /* general timer start */
  double elapsed_hash = 0;
  double elapsed_extract_dist = 0;
  double elapsed_mpi_send = 0;
  size_t nmsgs_sent = 0;
  u64 while_ctr = 0;
  char file_name[FILE_NAME_MAX_LENGTH];
  snprintf(file_name, sizeof(file_name), "data/stats/sender_%d", myrank);
  FILE* fp_timing = fopen(file_name, "a");

  /* this has to be a power of 2 please! */
  const size_t print_interval = (1LL<<25) - 1;
  
  u32* states_avx;
  int n_dist_points = 0; /* At the beginning we no dist points */
  size_t server_id=0, local_idx=0;

  /* clean the buffers using them! */
  memset(work_buf, 0, nbytes_per_server*NSERVERS); 
  memset(servers_counters, 0, sizeof(size_t)*NSERVERS);


  while (1) {
    /* increment the message counters after hashing */
    for (int lane = 0; lane<16; ++lane)
      ((u64*) Mavx[lane])[0] += 16;

    /* hash 16 messages */
    timer = wtime();
    states_avx = sha256_multiple_x16(Mavx);
    elapsed_hash += wtime() - timer;
    
    timer = wtime();
    extract_dist_points_16(states_avx,
			   Mavx,
			   digests,
			   msg_ctrs,
			   &n_dist_points);
    elapsed_extract_dist += wtime() - timer;
    
    
    /* put the distinguished points in specific serverss buffer */
    for (int i = 0; i<n_dist_points; ++i){ /* n_dist_points might be 0 */
      server_id = to_which_server(&digests[i*N]);

      /* where to put digest in server_id buffer? (this is a local view) */
      local_idx = servers_counters[server_id];

      /* 1 element = ctr||dgst  */
      /* #elments in  server buffer = PROCESS_QUOTA */
      /* There are NSERVERS in total */
      /* copy the ctr to the server buffer */
      memcpy(&work_buf[nbytes_per_server*server_id
		       + local_idx*one_elm_size],
	     &msg_ctrs[i], /* counter of the message i  */
	     sizeof(CTR_TYPE));

      /* copy the digest to the server buffer */
      memcpy(&work_buf[nbytes_per_server*server_id
		       + local_idx * one_elm_size
		       + sizeof(CTR_TYPE)], /* write after the counter part */
	     &digests[N*i],
	     N);


      /* we added an element in server n. server_id buffer */
      ++servers_counters[server_id];

      /* if the server buffer is full send immediately */
      if (servers_counters[server_id] == PROCESS_QUOTA){
	timer = wtime();
	buffered_isend(snd_bufs,
		       requests,
		       statuses,
		       nbufs,
		       &work_buf[server_id*nbytes_per_server],
		       nbytes_per_server,
		       server_id,
		       TAG_SND_DGST,
		       inter_comm);
	elapsed_mpi_send += wtime() - timer;
	++nmsgs_sent;
	
	servers_counters[server_id] = 0;
      }
    } /* end for*/
    ++while_ctr;
    
    if ( (while_ctr&print_interval) == 0){ /* while_ctr = 0 mod 2^25 */
      total_elapsed = wtime() - total_elapsed;
      
      show_and_save_benchmark(total_elapsed,
			      elapsed_mpi_send,
			      elapsed_hash,
			      elapsed_extract_dist,
			      nmsgs_sent,
			      nbytes_per_server,
			      16*while_ctr,
			      print_interval,
			      nsenders,
			      myrank,
			      fp_timing);
      
      
      /* reset all counters */
      total_elapsed = wtime();
      elapsed_extract_dist = 0;
      elapsed_hash = 0;
      elapsed_mpi_send = 0;
      nmsgs_sent = 0;
      while_ctr = 0;
    }
  } /* end while */
  
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
  u8 digests[16 * N] ; /* save the distinguished digests here (manually) */
  CTR_TYPE msg_ctrs[16]; /* save thde counters of the dist digests (manually) */

  

  /* Sending related buffers allocation  */

  /* |dgst| + |ctr| */
  size_t const one_pair_size = sizeof(u8)*N + sizeof(CTR_TYPE); 
  size_t const  nbytes_per_server = one_pair_size * PROCESS_QUOTA;

  /* pair = (ctr||dgst) */
  /* { (server0 paris..) | (server1 pairs...) | ... | (serverK pairs...) } */
  /* this where we put results */
  u8* work_buf = (u8*) malloc(nbytes_per_server * NSERVERS );

  
  /* this buffer will hold the message to be sent  */
  const int nbufs = 10;
  u8* snd_buf = (u8*) malloc(nbufs*nbytes_per_server);
  MPI_Status statuses[nbufs];
  MPI_Request* requests = (MPI_Request*) malloc(sizeof(MPI_Request)*nbufs);
  /* init requests */
  for (int i = 0; i<nbufs; ++i)
    requests[i] = MPI_REQUEST_NULL;

  
  /* ith_entry : how many non-sent element stored in server i buffer */
  size_t* server_counters = malloc(sizeof(size_t)*NSERVERS);


  if (snd_buf == NULL)
    puts("bsnd_buf is NULL");

  if (work_buf == NULL)
    puts("work_buf is NULL");


  /* memset(bsnd_buf, 0, (nbytes_per_server + MPI_BSEND_OVERHEAD)* NSERVERS); */
  memset(work_buf, 0, nbytes_per_server*NSERVERS);
  memset(server_counters, 0, sizeof(size_t)*NSERVERS);
  
  
  // ----------------------------- PART 1 --------------------------------- //
  // Regenrate the long message in parallel!                                //
  regenerate_long_message_digests(Mavx,
				  tr_states,
				  digests,
				  msg_ctrs,
				  work_buf,
				  snd_buf,
				  requests,
				  statuses,
				  nbufs,
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
  /* and set the initial counter */
  for (int lane = 0; lane<16; ++lane){
    memcpy(Mavx[lane], M, HASH_INPUT_SIZE);
    ((u64*) Mavx[lane])[0] += lane;
  }

  // print the template. this is not necessary.
  char txt[50];
  snprintf(txt, sizeof(txt), "sender #%d template=", myrank);
  print_byte_txt(txt, M,HASH_INPUT_SIZE);
  puts("\n");

  
  // ----------------------------- PART 2.b ---------------------------------- //
  // 1-  Sen the initial input to all receiving servers
  /* receivers before this stage are accepting any tag and and any source */
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Allgather(M,
		HASH_INPUT_SIZE,
		MPI_UNSIGNED_CHAR,
		NULL,
		0,
		MPI_UNSIGNED_CHAR,
		inter_comm);

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
			  snd_buf,
			  requests,
			  statuses,
			  nbufs,
			  server_counters,
			  nbytes_per_server,
			  one_pair_size,
			  myrank,
			  nsenders,
			  inter_comm);



  free(work_buf);
  free(snd_buf);
  free(server_counters);
  free(ctr_pt);
  

  return; // au revoir.

} // MPI_Finalize


