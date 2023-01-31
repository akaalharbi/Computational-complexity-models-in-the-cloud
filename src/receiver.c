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
  size_t ndigests = get_file_size(fp) / (N-DEFINED_BYTES);
  size_t nmemb = (ndigests/nchunks >=  1) ? ndigests/nchunks : 1;

  /* read digests from file to this buffer  */
  u8* digests = (u8*) malloc( (N-DEFINED_BYTES) * nmemb * sizeof(u8));


  /*  load one chunk each time */
  for (size_t i = 0; i<nchunks; ++i) {
    fread(digests, (N-DEFINED_BYTES), nmemb, fp);

    /* add them to dictionary */
    for (size_t j=0; j<nmemb; ++j) 
      dict_add_element_to(d,
			  &digests[j*(N-DEFINED_BYTES)] );
  }

  /* ndigests = nchunks*nmemb + remainder  */
  // read the remainder digests, since ndigests may not mutlipe of nchunks
  u8 stream_pt[N-DEFINED_BYTES];

  /* add as many hashes as possible */
  while ( !feof(fp) ){
    // use fread with a larger buffer @todo
    fread(stream_pt, sizeof(u8), (N-DEFINED_BYTES), fp);
    /* it adds the hash iff nprobes <= NPROBES_MAX */
    dict_add_element_to(d, stream_pt);
  }

  free(digests);
  return;
}


static inline void receiver_process_get_template(int myrank, int nproc, int nsenders, u8* templates)
{

  //---------------------------------------------------------------------------+
  // --- : receive the initial inputs from all generating processors 
  //---------------------------------------------------------------------------+

  // process : 0 -> NSERVERS - 1 (receivers),
  //         : NSERVERS -> nproc (senders)

  MPI_Status status;

  for (int i=0; i<nsenders; ++i){
    //printf("receiver %d is listening to %d for the template\n", myrank, i + NSERVERS);
    MPI_Recv(&templates[i*HASH_INPUT_SIZE],
	     HASH_INPUT_SIZE,
	     MPI_UNSIGNED_CHAR,
	     i + NSERVERS, /*  */
	     TAG_RANDOM_MESSAGE,
	     MPI_COMM_WORLD,
	     &status);

  }
}





void receiver_process_task(dict* d,
			   int myrank,
			   int nproc,
			   int nsenders,
			   u8* templates )
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
  FILE* fp = fopen(file_name, "a");

  int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                    + sizeof(CTR_TYPE); /* |dgst| + |ctr| */
  
  //+ receive message in this buffer
  u8* rcv_buf = (u8*) malloc(one_pair_size * PROCESS_QUOTA);
  
  //+ copy rcv_buf to lookup_buf then listen to other sender using rcv_buf
  u8* lookup_buf = (u8*) malloc(one_pair_size * PROCESS_QUOTA);

  size_t rcv_array_size = one_pair_size*PROCESS_QUOTA;
  /* u8* initial_inputs = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nproc_snd); */
  
  MPI_Status status;
  MPI_Request request;

  int* indices = (int*)malloc(sizeof(int)*nproc);

  /* How many candidates were stored? and remove partial candidates */
  size_t nfound_cnd = get_file_size(fp) / HASH_INPUT_SIZE ;
  // truncate the candidates file, in case a partial candidate was stored
  truncate(file_name, nfound_cnd*HASH_INPUT_SIZE);
  size_t old_nfound_candidates = nfound_cnd;

  /* 1st sender has rank = NSERVER -scaling-> 1st sender name = 0 */
  int sender_name_scaled = 0; 


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
	   MPI_COMM_WORLD,
	   &status);


  printf("recv #%d got its 1st message from %d\n",
	 myrank,
	 status.MPI_SOURCE);
  
  // copy the received message and listen immediately
  memcpy(lookup_buf, rcv_buf, rcv_array_size);
  /* print_char(lookup_buf, rcv_array_size); */
  /* puts("-=-=-=-=-=-=-=-=-=-=-=-=-="); */
  
  /* 1st sender has rank = NSERVER -scaling-> 1st sender name = 0 */
  sender_name_scaled = status.MPI_SOURCE - NSERVERS; // who sent the message?


  while (NNEEDED_CND > nfound_cnd) {    
    /* printf("nfound_cnd = %lu, myrank=%d\n", nfound_cnd, myrank); */
    //+ receive messages from different processors
    MPI_Irecv(rcv_buf, /* store in this location */
	      rcv_array_size, 
	      MPI_UNSIGNED_CHAR,
	      MPI_ANY_SOURCE, /* sender */
	      TAG_SND_DGST, 
	      MPI_COMM_WORLD,
	      &request);


    // probe these messages and update the founded candidates
    /* printf("recv#%d probing sender #%d messages\n", */
    /* 	   myrank, sender_name_scaled); */
    nfound_cnd += lookup_multi_save(d, /* dictionary to look inside */
				    lookup_buf, /* messages to search in d */
				    &templates[sender_name_scaled
					       *HASH_INPUT_SIZE],
				    PROCESS_QUOTA,/* how many msgs in rcv_buf */
				    fp, /* file to record cadidates */
				    myrank,
				    sender_name_scaled + NSERVERS);


    if (nfound_cnd - old_nfound_candidates > 0) {
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	     "receiver #%d has %lu out of %llu candidates from #%d\n"
	     "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n",
	     myrank, nfound_cnd, NNEEDED_CND, status.MPI_SOURCE);
      
      old_nfound_candidates = nfound_cnd;
    }
    MPI_Wait(&request, &status);


   
    sender_name_scaled = status.MPI_SOURCE - NSERVERS; // get the name of the new sender
    memcpy(lookup_buf, rcv_buf, rcv_array_size);

    /* printf("from %d:\n", status.MPI_SOURCE); */
    /* print_char(lookup_buf, rcv_array_size); */
    /* print_attack_information(); */

    }

  // good job
  free(rcv_buf);
  free(indices);
  fclose(fp);

  printf("recv #%d done a good job\n", myrank);

  exit(EXIT_SUCCESS);
}


void receiver(int myrank, MPI_Comm mpi_communicator, int nsenders)
{
  //--------------------------------------------------------------------------+
  // I'm a receiving process: receive hashes, probe them, and send candidates |
  // Process Numbers: [0,  NSERVERS - 1]                                      |
  //--------------------------------------------------------------------------+

  int nproc;
  MPI_Comm_size(mpi_communicator, &nproc);  
  u8* templates = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nsenders);

  // PART 1: LOAD dictionary
  char txt[100];
  snprintf(txt, sizeof(txt), "recv#%d before dict load", myrank);
  print_memory_usage(txt);

  dict* d = dict_new(NSLOTS_MY_NODE);
  double time_start = wtime();
  char file_name[FILE_NAME_MAX_LENGTH];
  snprintf(file_name, sizeof(file_name), "data/digests/%d", myrank);
  FILE* fp = fopen(file_name, "r");
  load_file_to_dict(d, fp);


  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	 "recv #%d dict read in %0.2fsec\n"
	 "It has %lu elms, file has %lu elms\n"
	 "d->nslots = %lu, d->nelements=%lu, filling rate=%f \n"
	 "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
	 myrank, wtime() - time_start,
	 d->nelements, d->nelements_asked_to_be_inserted,
	 d->nslots, d->nelements,
	 ((float) d->nelements)/d->nslots);

  snprintf(txt, sizeof(txt), "recv#%d after  dict load", myrank);
  print_memory_usage(txt);
  fclose(fp);


  
  // PART 2: Get templates from all senders
  receiver_process_get_template(myrank, nproc, nsenders, templates);

  // PART 3: Get templates from all senders
  /* listen to senders, probe their digest, in case a candidate is found: */
  /* save the message the generates the cadidate digest. */
  /* exit when enough number of candidates are found i.e. NNEEDED_CND */
  receiver_process_task(d, myrank, nproc, nsenders, templates);


  free(templates);
}
