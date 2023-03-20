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
#include <unistd.h>
#include "common.h"
#include "sender.h"




//---------------------------- UTILITY FUNCTIONS -------------------------------
// local function


static inline int lookup_multi_save(dict *d,
				    u8 *ctrs_dgsts, /* known bytes are removed  */ 
                                    u8 template[HASH_INPUT_SIZE], /* initial random message */
                                    size_t npairs, /* #pairs in rcv_buf */
				    FILE *fp, /* store candidates here */
				    int myrank, /* for debugging */
				    int source) /* for debugging */
{
  // -------------------------------------------------------------------------+
  // Given ctrs_dgsts={ctr1||dgst1,..} from specific rank, probes each dgst,if|
  // prope(dgst)=/=0 store its related msg in fp, and returns the number of   |
  // stored messages.                                                         |
  // -------------------------------------------------------------------------+
  // INPUT:                                                                   |
  // - dict                                                                   |
  // - templaet: random message with zero counter                             |
  // - *ctrs_dgsts: array of pairs (ctr, digest ) packed as ctr||digest       |
  // - npairs: how many pairs (ctr, dgst) in the stream.                      |
  // - fp: where to store candidate messages, i.e. messages whose digest      |
  //       returns positive answser when it gets probed in the dictionary     |
  // - myrank: myrank as an MPI, process (we don't use it!)                   |
  // - source: who sent the message                                           |
  // -------------------------------------------------------------------------+
  // NOTES:                                                                   |
  // |msg| = sizeof(CTR_TYPE) since we receive the counter                    |
  // |full-msg| = NWORDS_INPUT*WORD_SIZE (config.h) it can be constructed     |
  // using msg and template.                                                  |
  // |dgst| = N-DEFINED_BYTES (should be N,but we skip know bits .e.g nserver)|
  // - This function doesn't need the sender rank, it alread has the sender's |
  //   template.                                                              |
  // -------------------------------------------------------------------------+


  /* the stream size is multiple of one_pair size */
  static int one_pair_size = sizeof(CTR_TYPE) + (N-DEFINED_BYTES)*sizeof(u8);


  
  /* how many digests give positive answer they get probed  */
  int npositive_dgsts = 0;
  /* did dict say it has the digest? (it doesn't have to be correct)*/
  int is_positive = 0; 

  /* construct the full message here */
  static u8 M[HASH_INPUT_SIZE]; 
  
  for (size_t i=0; i<npairs; ++i){
    /* dictionary only read |dgst| bytes by default  */
    /* probing pair #i */

    /* pair := ctr||dgst, probe d for dgst  */
    is_positive =  dict_has_elm(d, /* it probes the digest */
			&ctrs_dgsts[i*one_pair_size /* Go 2 the ith pair */
				    + sizeof(CTR_TYPE)] /* skip ctr part */);
    
    
    if (is_positive){ /* did we find a candidate msg||dgst? */
      /* reconstructing the message: get the template */
      memcpy(M, template, HASH_INPUT_SIZE);

      /* set the counter part */
      memcpy(M,
	     &ctrs_dgsts[i*one_pair_size],
	     sizeof(CTR_TYPE));
      
      /* assert(is_dist_msg(M)); // debugging  */

      /* finally write the reconstructed message */
      fwrite(M, sizeof(u8), HASH_INPUT_SIZE, fp);
      fflush(fp); /* ensures it's written */
      ++npositive_dgsts;
    }
  }
  return npositive_dgsts;
}




// -----------------------------------------------------------------------------
// local function
// this function should be rewritten @todo 
static void write_digest_to_dict(dict *d,
				  u8* restrict rcv_buf,
				  const int nsenders,
				  MPI_Comm inter_comm)
{
  // ==========================================================================+
  // Summary: Listen to senders who regenerates the long messsage digests, and |
  //          record the received digests in the dictionary.                   |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d`: Dictionary that will keep elements from *fp.                        |
  // `rcv_buf`: this buffer will be used to record message from senders.       |
  // `add_buf`: this is used to copy the rcv_buf, so it can listen again.      |
  // `nsenders`: this is used to know when to stop listening. since each sender|
  //             after done hashing, it will send a different tag.             |
  // `inter_comm`: MPI inter communication object that seperates senders and   |
  //               receivers into two disjoint groups.                         |
  // ==========================================================================+

  // @todo note: not to forget we will change the direction right most bytes
  // contain the distinguished point zeros and the bits.
  MPI_Status status;
  /* A sender will send a 0 tag if it needs to hash more, otherwise it will send tag = 1 */
  int ncompleted_senders = 0; 
  size_t rcv_size = PROCESS_QUOTA*N;


  //---------------------------------------------------------------------------+
  // Receive digests and add them to dictionary
  
  /* Receive one message before */
  while (ncompleted_senders < nsenders) {
    MPI_Recv(rcv_buf, /* store in this location */
	      rcv_size, /* How many bytes to receive  */
	      MPI_UNSIGNED_CHAR,
	      MPI_ANY_SOURCE, /* any sender */
	      MPI_ANY_TAG,  /* tag = 1 means a sender has done its work */ 
	      inter_comm,
	      &status);


    /* add them to dictionary */
    for (size_t j=0; j<PROCESS_QUOTA; ++j) 
      dict_add_element_to(d, &rcv_buf[N*j]);

    /* If a sender is done hashing, it will make rcv_buf = {0}, and has tag=1 */
    /* The dictionary by design will ignore all zero messages */
    ncompleted_senders += status.MPI_TAG;
  }
}





void receiver_process_task(int const myrank,
			   size_t const rcv_array_size,
			   dict* restrict d,
			   u8* restrict templates,
			   u8* restrict rcv_buf,
			   u8* restrict lookup_buf,
			   MPI_Comm inter_comm)
{
  // todo check the loops, currently they are errornous!
  //---------------------------------------------------------------------------+
  // ------- I am a receiving processor, I only probe the dictionary      -----|
  // Receive messages till we found at least NNEEDED_CND_THIS_SERVER messages  |
  // candidates. Send candidates to archive process.                           |
  //---------------------------------------------------------------------------+
  // nproc_snd: number of sender processes
  // nproc : number of all processes

  //---------------------------------------------------------------------------+
  // ------- Part 1 : Init receiving buffers of digests and counters      -----|
  //---------------------------------------------------------------------------+
  /* create file: data/messages/myrank that will hold messages whose hashes */
  /* gives a postivie response when probing the dictionary */
  
  char file_name[FILE_NAME_MAX_LENGTH]; /* "data/messages/%d" */
  snprintf(file_name, sizeof(file_name), "data/messages/%d", myrank );
  FILE* fp = fopen(file_name, "a"); /* register message candidates here */
  
  MPI_Status status;
  MPI_Request request;
  /* MPI_Request_free(&request); /\* in 1st time there is no waiting  *\/ */
  
  /* How many candidates were stored? and remove partial candidates */
  size_t nfound_cnd = get_file_size(fp) / HASH_INPUT_SIZE ;
  // truncate the candidates file, in case a partial candidate was stored
  truncate(file_name, nfound_cnd*HASH_INPUT_SIZE);
  size_t old_nfound_candidates = nfound_cnd;

  /* with inter-communication, 1st sendre has name 0 */
  int sender_name = 0;  


  //---------------------------------------------------------------------------+
  // --- Part 3: Receive digests and probe them
  //---------------------------------------------------------------------------+
  /* printf("recv#%d is listening\n", myrank); */
  // listen the first time

  printf("recv #%d will be posted \n", myrank);
  MPI_Recv(rcv_buf,
	   rcv_array_size,
	   MPI_UNSIGNED_CHAR,
	   MPI_ANY_SOURCE,
	   TAG_SND_DGST,
	   inter_comm,
	   &status);


  printf("recv #%d got its 1st message from %d\n",
	 myrank,
	 status.MPI_SOURCE);
  
  // copy the received message and listen immediately
  memcpy(lookup_buf, rcv_buf, rcv_array_size);
  /* print_char(lookup_buf, rcv_array_size); */
  /* puts("-=-=-=-=-=-=-=-=-=-=-=-=-="); */
  
  /* The sender rank in its respective local group  */
  sender_name = status.MPI_SOURCE; // - NSERVERS; // who sent the message?


  while (NNEEDED_CND > nfound_cnd) {    
    /* printf("nfound_cnd = %lu, myrank=%d\n", nfound_cnd, myrank); */
    //+ receive messages from different processors
    MPI_Irecv(rcv_buf, /* store in this location */
	      rcv_array_size, 
	      MPI_UNSIGNED_CHAR,
	      MPI_ANY_SOURCE, /* sender */
	      TAG_SND_DGST, 
	      inter_comm,
	      &request);

    /* printf("recv#%d probing sender #%d messages\n", */
    /* 	   myrank, sender_name_scaled); */

    /* probe these messages and update the founded candidates */
    nfound_cnd += lookup_multi_save(d, /* dictionary to look inside */
				    lookup_buf, /* messages to search in d */
				    &templates[sender_name * HASH_INPUT_SIZE],
				    PROCESS_QUOTA,/* how many msgs in rcv_buf*/
				    fp, /* file to record cadidates */
				    myrank, /* was-this only for debugging? */ // @todo
				    sender_name); /* why do we need sender name here? */


    if (nfound_cnd - old_nfound_candidates > 0) {
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	     "receiver #%d has %lu out of %llu candidates from #%d\n"
	     "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n",
	     myrank, nfound_cnd, NNEEDED_CND, status.MPI_SOURCE);
      
      old_nfound_candidates = nfound_cnd;
    }
    MPI_Wait(&request, &status);


    /* update buffers according to the new message */
    sender_name = status.MPI_SOURCE; // get the name of the new sender
    memcpy(lookup_buf, rcv_buf, rcv_array_size);

    /* printf("from %d:\n", status.MPI_SOURCE); */
    /* print_char(lookup_buf, rcv_array_size); */
    /* print_attack_information(); */

    }

  // good job
  free(rcv_buf);
  /* free(indices); */
  fclose(fp);

  printf("recv #%d done a good job\n", myrank);

  exit(EXIT_SUCCESS);
}


void receiver(int local_rank, 
	      int nsenders, /* How many senders in the remote group  */
	      MPI_Comm inter_comm)
{
  //--------------------------------------------------------------------------+
  // I'm a receiving process: receive hashes, probe them, and send candidates |
  // Process Numbers: [0,  NSERVERS - 1]                                      |
  // inter_comm: for communicating with senders/                              |
  // locaL_comm: essentially there is no communications going in this comm    |
  //--------------------------------------------------------------------------+


  //------------------------ INIT VARIABLES ----------------------------------+

  int const  one_pair_size = sizeof(u8)*N
                           + sizeof(CTR_TYPE); /* |dgst| + |ctr| */

  size_t const rcv_array_size = one_pair_size*PROCESS_QUOTA;
  
  // receive message in this buffer regardless if it's pure digests or
  // ctr||digests 
  u8* rcv_buf = (u8*) malloc(one_pair_size * PROCESS_QUOTA);
  //+ copy rcv_buf to lookup_buf then listen to other sender using rcv_buf
  u8* lookup_buf = (u8*) malloc(one_pair_size * PROCESS_QUOTA);

  memset(rcv_buf, 0, one_pair_size*PROCESS_QUOTA);
  memset(lookup_buf, 0, one_pair_size*PROCESS_QUOTA);


  /* save initial random messages of each sender here  */
  u8* templates = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nsenders);


  char txt[100];
  snprintf(txt, sizeof(txt), "recv#%d before dict load", local_rank);
  print_memory_usage(txt);

  dict* d = dict_new(NSLOTS_MY_NODE);
  double time_start = wtime();


  //--------------------------------- PART 1 ----------------------------------+
  // PART 1: Get the long message digests from all senders
  write_digest_to_dict(d,
		       rcv_buf,
		       nsenders,
		       inter_comm);


  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	 "recv #%d dict read in %0.2fsec\n"
	 "It has %lu elms, file has %lu elms\n"
	 "d->nslots = %lu, d->nelements=%lu, filling rate=%f \n"
	 "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
	 local_rank, wtime() - time_start,
	 d->nelements, d->nelements_asked_to_be_inserted,
	 d->nslots, d->nelements,
	 ((float) d->nelements)/d->nslots);

  snprintf(txt, sizeof(txt), "recv#%d after  dict load", local_rank);
  print_memory_usage(txt);



  //--------------------------------- PART 2 ----------------------------------+
  // corresponds to part 2.b in senders
  // PART 2: Get templates from all senders
  MPI_Allgather(NULL, 0, NULL, /* receivers don't send */
		templates, /* save messages here */
		HASH_INPUT_SIZE,
		MPI_UNSIGNED_CHAR,
		inter_comm);


  //--------------------------------- PART 2 ----------------------------------+
  /* listen to senders, probe their digest, in case a candidate is found: */
  /* save the message the generates the cadidate digest. */
  /* exit when enough number of candidates are found i.e. NNEEDED_CND */
  receiver_process_task(local_rank,
			rcv_array_size,
			d, /* dictionary pointer  */
			templates,
			rcv_buf,
			lookup_buf,
			inter_comm);


  free(templates);
}
