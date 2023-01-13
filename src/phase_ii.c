// phase ii: high level overview 
// four types of processors: senders (the majority), receivers (#NSERVERS),
// Upon starting the program:
//
// `receiver`:- (init) load digests into a dictionary(once), receive messages
//              templates from senders (once).
//            - (listen)recieve digests from senders, probe them, and store the
//               candidates.
//            - (send to archive) if a receiver had enough digests, send them to
//              to the archive. Then go back to the listen state.
//
// `sender`:- (init) generate random message templates, send it to all recivers
//          - (gen) hash many random messages, decides which receiver x should
//              get a sepcific digest, store the digest in a buffer x, if it
//              is reach the quota, send that buffer to receiver x.
//              repeat infinitely
//





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


//
// checklist: msg||dgst 

//---------------------------- UTILITY FUNCTIONS -------------------------------






static inline int lookup_multi_save(dict *d,
				    u8 *ctrs_dgsts, /* known bytes are removed  */ 
                                    u8 template[HASH_INPUT_SIZE], /* initial random message */
                                    size_t npairs, /* #pairs in rcv_buf */
				    FILE *fp, /* store candidates here */
				    int myrank, /* for debugging */
				    int source) /* for debugging */
{
  // -------------------------------------------------------------------------+
  // Given ctrs_dgsts={ctr1||dgst1,..} from specific rank, probes each dgst, if   |
  // prope(dgst)=/=0 store its related msg in fp, and returns the number of   |
  // stored messages.                                                         |
  // -------------------------------------------------------------------------+
  // INPUT:                                                                   |
  // - dict                                                                   |
  // - stream                                                                 |
  // - npairs: how many pairs (ctr, dgst) in the stream.                      |
  // -------------------------------------------------------------------------+
  // NOTES:                                                                   |
  // |msg| = sizeof(CTR_TYPE) since we receive the counter                    |
  // |full-msg| = NWORDS_INPUT*WORD_SIZE (config.h) it can be constructed     |
  // using msg and template.                                                  |
  // |dgst| = N-DEFINED_BYTES (should be N,but we skip know bits .e.g nserver)|
  // -------------------------------------------------------------------------+
  /* the stream size is multiple of one_pair size */
  static int one_pair_size = sizeof(CTR_TYPE) + (N-DEFINED_BYTES)*sizeof(u8);
  static int msg_size = HASH_INPUT_SIZE; /* full message not only the counter */

  
  /* how many digests give positive answer they get probed  */
  int npositive_dgsts = 0; 
  int is_positive = 0; /* what is this?  */

  /* construct the full message here */
  static u8 M[HASH_INPUT_SIZE]; 
  
  for (size_t i=0; i<npairs; ++i){
    /* dictionary only read |dgst| bytes by default  */
    /* probing pair #i */

    is_positive =  dict_has_elm(d,
			&ctrs_dgsts[i*one_pair_size /* Go 2 the ith pair */
				    + sizeof(CTR_TYPE)] /* skip ctr part */
				);
    

    
    if (is_positive){ /* did we find a candidate msg||dgst? */

      /* reconstructing the message: get the template */
      memcpy(M, template, HASH_INPUT_SIZE);

      /* set the counter part */
      memcpy(M,
	     &ctrs_dgsts[i*one_pair_size],
	     sizeof(CTR_TYPE));

      /* assert(is_dist_msg(M)); */ // debugging 

      /* finally write the reconstructed message */
      fwrite(M, sizeof(u8), msg_size, fp);
      fflush(fp); /* ensures it's written */
      ++npositive_dgsts;
    }
  }
  return npositive_dgsts;
}


// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
void load_file_to_dict(dict *d, FILE *fp)
{
  // ==========================================================================+
  // Summary: Load hashes from the file *fp, and try to store them in dict *d  |
  // Note: dict d has the right to reject inserting element.                   |
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

  
  size_t nchunks = 1000; 
  size_t ndigests = get_file_size(fp) / (N-DEFINED_BYTES);
  size_t nmemb = (ndigests/nchunks >=  1) ? ndigests/nchunks : 1;
  size_t remainder = ndigests % nchunks;


  u8* digests = (u8*) malloc( (N-DEFINED_BYTES) * nmemb * sizeof(u8));


  if(nmemb != (ndigests/nchunks)){
    // the file has less than 1000 digests
    u8 stream_pt[N-DEFINED_BYTES];
    /* add as many hashes as possible */
    while ( !feof(fp) ){
      // use fread with a larger buffer @todo
      fread(stream_pt, sizeof(u8), (N-DEFINED_BYTES), fp);
      /* it adds the hash iff nprobes <= NPROBES_MAX */
      dict_add_element_to(d, stream_pt);
    }
  } 
  
  for (size_t i = 0; i<nchunks; ++i) {
    fread(digests, (N-DEFINED_BYTES), nmemb, fp);

    /* add them to dictionary */
    for (size_t j=0; j<nmemb; ++j) 
      dict_add_element_to(d,
			  &digests[j*(N-DEFINED_BYTES)] );
    
  }

  // read the left digests 
  fread(digests, (N-DEFINED_BYTES), remainder, fp);

  /* add them to dictionary */
  for (size_t j=0; j<remainder; ++j) 
    dict_add_element_to(d,
			&digests[j*(N-DEFINED_BYTES)] );


  
  // fclose(fp); // don't close the file 
}


void send_random_message_template(u8 M[HASH_INPUT_SIZE])
{ /* Send M immediately and clear the memory at the end  */
  
  MPI_Request requests[NSERVERS]; /* they will be used only here  */
  MPI_Status statuses[NSERVERS];
  // how about collective communications? I can't find a simple way using them.
  for (int i=0; i<NSERVERS; ++i) {
    MPI_Isend(M,/* snd_buf */
	      HASH_INPUT_SIZE,
	      MPI_UNSIGNED_CHAR,
	      i /* to whom */,
	      TAG_RANDOM_MESSAGE,
	      MPI_COMM_WORLD, &requests[i]);
  }
  MPI_Waitall(NSERVERS, requests, statuses);
} /* clear stack variables */



void sender(int myrank, MPI_Comm mpi_communicator)
{

  // ------------------------------------------------------------------------+
  // I am a sending processor. I only generate hashes and send them.
  // Process Numbers: [NSERVERS + 1,  nproc]
  //-------------------------------------------------------------------------+

  /* init state: this should goes inside sender_process_task */
  /* set up a random message */

  /* M = 64bit ctr || 64bit nonce || random value */
  u8 M[HASH_INPUT_SIZE]; /* random word */

  /* Get a random message only once */
  CTR_TYPE* ctr_pt = (CTR_TYPE*) M; /* counter pointer  */
  getrandom(M, HASH_INPUT_SIZE, 1);
  ctr_pt[0] = 0; /* zeroing the first 64bits of M */

  char txt[50];
  snprintf(txt, sizeof(txt), "sender #%d template=", myrank);
  print_byte_txt(txt, M,HASH_INPUT_SIZE);
    
  /* Send the initial input to all receiving servers */
  /* printf("sender #%d is going to send its template\n", myrank); */
  send_random_message_template(M); 
  /* printf("sender #%d done sending template!\n", myrank); */
  /* printf("sender #%d template=\n", myrank); */
  /* print_char(M, 64); */

  /* generate hashes and send them to a servers */


  /* generate hashes and send them to a server as soon its buffer is complete */
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}
  // dgst := N bytes of the digest
  // ctr  := (this is a shortcut that allows us not to send the whole message)
  //---------------------------------------------------------------------------+
  WORD_TYPE Mstate[NWORDS_STATE] = {HASH_INIT_STATE};
  size_t one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                    + sizeof(CTR_TYPE); /* |dgst| + |ctr| - |known bits|*/

  CTR_TYPE msg_ctr = 0;
  int server_number;


  // ---- Part1 : initializing sending buffers 
  /* int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES) */
  /*                      + sizeof(CTR_TYPE); /\* |dgst| + |ctr| *\/ */

  u8* snd_buf = (u8*) malloc(one_pair_size
			     *PROCESS_QUOTA
			     *NSERVERS);


  /* Decide where to place the kth digest in server i buffer */
  // How many bytes are reserverd for each server in snd_buf
  size_t nbytes_per_server = one_pair_size * PROCESS_QUOTA;
  size_t offset = 0; /* where inside the buffer should we go */
  const u64 mask_test = (1LL<<DIFFICULTY) - 1;

  /* pos i: How many messages we've generated to be sent to server i? */
  u64 servers_ctr[NSERVERS] = {0};
  int snd_ctr = 0; /* how many messages have been sent */
  double time_end = wtime();
  /* set number of generated messages before sending to 0 */
  /* memset(servers_ctr, 0, sizeof(u32)*NSERVERS);  */

  //-------------------------------------------------------------------------+
  //------- Part 2 : Generate hashes and send them  -------
  //-------------------------------------------------------------------------+
 
  printf("sender #%d: done init mpi, and sharing the its template."
	 "going to generate hashes \n", myrank);
  
  // find_hash_distinguished_init(); 



  while(1) { /* when do we break? never! */
  /* while(i<1) { /\* when do we break? never! *\/ */

    /* Find a message that produces distinguished point */
    find_hash_distinguished( M, Mstate, &msg_ctr, mask_test);

    //+ decide to which server to add to? 
    server_number = to_which_server((u8*) Mstate);


    /* 1st term: go to server booked memory, 2nd: location of 1st free place*/
    offset = server_number * nbytes_per_server
           + servers_ctr[server_number] * one_pair_size;
    // recall que one_pair_size =  |dgst| + |ctr| - |known bits|

    /* printf("sender #%d, server_ctr=0x%llu, server_number=%d\n", */
    /* 	   myrank, servers_ctr[server_number], server_number); */

    // record a pair (msg, dgst), msg is just the counter in our case
    /* record the counter  */
    
    memcpy(&snd_buf[ offset ],
	   M,
	   sizeof(CTR_TYPE) );

    /* After the counter save N-DEFINED_BYTES of MState */
    memcpy( &snd_buf[offset + sizeof(CTR_TYPE)], /* copy digest to snd_buf[offset] */
	    ((u8*)Mstate) + DEFINED_BYTES, /* skip defined bytes */
	    N-DEFINED_BYTES );

    servers_ctr[server_number] += 1;
    

    /* if (server_number == 2){ */
      

    /*   printf("server=%d, server_ctr=%llu,\n",server_number, */
    /* 	     servers_ctr[server_number]); */
    /*   print_char(&snd_buf[offset], one_pair_size); */


    /*   /\* this server has one more digest *\/ */

    /*   for (int i=0; i<NSERVERS; ++i) { */
    /* 	printf("%d_ctr=%llu, offset=%lu, nbytes/server=%lu, one_pair %lu\n", */
    /* 	       i, servers_ctr[i], offset, nbytes_per_server, one_pair_size); */
    /*   } puts(""); */
    /* } */
    
    if (servers_ctr[server_number] >= PROCESS_QUOTA){
      
      /* printf("rank %d: sending to %d, snd_ctr=%d, it took %fsec\n", */
      /* 	     myrank, server_number, snd_ctr, wtime() - time_start); */
      /* time_start = wtime(); */
      /* we have enough messages to send to server (server_number) */

      /* printf("===============================================\n" */
      /*        "sender #%d -> recv #%d before sending %0.2fsec\n" */
      /* 	     "===============================================\n\n", */
      /* 	     myrank, server_number, wtime() -  time_end ); */

      MPI_Send(&snd_buf[server_number*nbytes_per_server],
		PROCESS_QUOTA*one_pair_size,
		MPI_UNSIGNED_CHAR,
		server_number,
		TAG_SND_DGST,
		MPI_COMM_WORLD);

      ++snd_ctr;

      /* printf("-----------------------------------------------\n" */
      /*        "sender #%d -> recv #%d sending done %0.2fsec\n" */
      /* 	     "-----------------------------------------------\n\n", */
      /* 	     myrank, server_number, wtime() -  time_end ); */

      /* 	snprintf(txt, sizeof(txt), "sender#%d,  server=%d, snd_buf=", */
      /* 		 myrank, server_number ); */
      /* 	print_byte_txt(txt, */
      /* 		       &snd_buf[server_number*nbytes_per_server], */
      /* 		       one_pair_size*PROCESS_QUOTA); */
      /* printf("rank %d: sending done to %d, snd_ctr=%d, it took %fsec\n", */
      /* 	     myrank, server_number, snd_ctr, wtime() - time_start); */
      /* time_start = wtime(); */
      // todo do we need here buffer detach? 
      /* It is enough to reset the counter. The memroy will be rewritten */  
      servers_ctr[server_number] = 0;
      time_end = wtime();
    }
  } 


  free(snd_buf);
  return; // au revoir.

} // MPI_Finalize



static inline void receiver_process_get_template(int myrank, int nproc, int nproc_snd, u8* templates)
{



  //---------------------------------------------------------------------------+
  // --- : receive the initial inputs from all generating processors 
  //---------------------------------------------------------------------------+

  // process : 0 -> NSERVERS - 1 (receivers),
  //         : NSERVERS -> nproc (senders)

  MPI_Status status;

  for (int i=0; i<nproc_snd; ++i){
    //printf("receiver %d is listening to %d for the template\n", myrank, i + NSERVERS);
    MPI_Recv(&templates[i*HASH_INPUT_SIZE],
	     HASH_INPUT_SIZE,
	     MPI_UNSIGNED_CHAR,
	     i + NSERVERS,
	     TAG_RANDOM_MESSAGE,
	     MPI_COMM_WORLD,
	     &status);

    //print_char(&templates[i*HASH_INPUT_SIZE], HASH_INPUT_SIZE);
    
    /* print_char(&initial_inputs[i*HASH_INPUT_SIZE], HASH_INPUT_SIZE); */
  }
}







void receiver_process_task(dict* d, int myrank, int nproc, int nproc_snd, u8* templates )
{
  // todo check the loops, currently they are errornous!
  //---------------------------------------------------------------------------+
  // ------- I am a receiving processor, I only probe the dictionary      -----|
  // Receive messages till we found at least NNEEDED_CND_THIS_SERVER messages  |
  // candidates. Send candidates to archive process.                           |
  //---------------------------------------------------------------------------+
  // nproc_snd: number of sender processes
  // nproc : number of all processes 
  // ------- Part 1 : Init receiving buffers of digests and counters      -----|
  //---------------------------------------------------------------------------+

  /* create file: data/messages/myrank that will hold messages whose hashes */
  /* gives a postivie response when probing the dictionary */
  
  char file_name[FILE_NAME_MAX_LENGTH]; /* "data/send/messages/%d" */
  snprintf(file_name, sizeof(file_name), "data/send/messages/%d", myrank );
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

  size_t nfound_cnd = get_file_size(fp) / HASH_INPUT_SIZE ;
  size_t old_nfound_candidates = nfound_cnd;
  /* 1st sender has rank = NSERVER -scaling-> 1st sender name = 0 */
  int sender_name_scaled = 0; 


  
  
  long vmrss_kb, vmsize_kb;
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  /* printf("reciver #%d memory usage after mpi init: ram: %ld kb vm: %ld kb\n", */
  /* 	 myrank, vmrss_kb, vmsize_kb); */
  




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

  printf("recv #%d got a message from %d\n"
	 "(We will message you again)\n", myrank, status.MPI_SOURCE);
  
  // copy the received message and listen immediately
  memcpy(lookup_buf, rcv_buf, rcv_array_size);
  /* 1st sender has rank = NSERVER -scaling-> 1st sender name = 0 */
  sender_name_scaled = status.MPI_SOURCE - NSERVERS; // who sent the message?


  // @todo restore the while loop
  /* while (NNEEDED_CND  > nfound_cnd) {  */ // @todo fix this line
  // strange behavior nneeded_cnd = 0
    /* from = status.MPI_SOURCE; */
    /* char txt[50]; */
    /* snprintf(txt, sizeof(txt), "rcv #%d, receive buffer=", myrank); */
    /* print_byte_txt(txt,rcv_buf, rcv_array_size); */





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
    
    //+ probe these messages and update the founded candidates
    /* printf("recv#%d probing sender #%d messages\n", */
    /* 	   myrank, sender_name_scaled); */
    nfound_cnd += lookup_multi_save(d, /* dictionary to look inside */
				    lookup_buf, /* messages to search in d */
				    &templates[sender_name_scaled
					       *HASH_INPUT_SIZE],
				    PROCESS_QUOTA,/* how many msgs in rcv_buf */
				    fp,
				    myrank,
				    status.MPI_SOURCE); /* file to record cadidates */

    if (nfound_cnd - old_nfound_candidates > 0) {
      printf("+++++++++++++++++++++++++++++++++++++++++\n"
	     "receiver #%d has %lu candidates from #%d\n"
	     "++++++++++++++++++++++++++++++++++++++++++\n\n",
	     myrank, nfound_cnd, status.MPI_SOURCE);
      
      old_nfound_candidates = nfound_cnd;
    }
    MPI_Wait(&request, &status);


      /* from = status.MPI_SOURCE; */
      /* char txt[50]; */
      /* snprintf(txt, sizeof(txt), "rcv #%d, from #%d, receive buffer=", */
      /* 	       from, */
      /* 	       myrank); */
      /* print_byte_txt(txt,rcv_buf, rcv_array_size); */
   
    sender_name_scaled = status.MPI_SOURCE - NSERVERS; // get the name of the new sender
    memcpy(lookup_buf, rcv_buf, rcv_array_size);
    }


  // good job
  free(rcv_buf);
  free(indices);
  fclose(fp);

  printf("recv #%d done a good job\n", myrank);
  /* exit(EXIT_SUCCESS); */
}



void receiver(int myrank, MPI_Comm mpi_communicator){
  
}


// 
// -----------------------------------------------------------------------------


int main(int argc, char* argv[])
{

  // the program is a bit stupid, it can't handle when N <= DEFINED_BYTES 
  assert(N - DEFINED_BYTES > 0);

  // some extra arguments to communicate with other
  
  // ==========================================================================+
  // Summary: Three main characters:                                           |
  //  1- producers: rank [NSERVERS+1, INF] generate hashes indifinitely.       |
  //  2- consumers: ranks [0, NSERVERS-1] receive hashes and probe them in the |
  //                in its dictionary. Save those that return positive answer  |
  //  3- archive: rank=NSERVER store all messages||hashes that dict(hash) > 0  |
  // ==========================================================================+



  // --------------------- INIT MPI & Shared Variables ------------------------+
  int nproc, myrank;
  
  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* How many procs that are going t send */
  int nproc_snd = nproc - NSERVERS;

  printf("There are %d senders\n", nproc_snd);
  // Who am I? a sender,  or a receiver?
  if (myrank >= NSERVERS){
    sender(myrank, MPI_COMM_WORLD); /* never ends :) */
  }
  else if (myrank < NSERVERS){ /* receiver, repeat infinitely  */

    long vmrss_kb, vmsize_kb;
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    printf("reciver #%d memory usage beginning: ram: %ld kb vm: %ld kb\n"
	   ,myrank, vmrss_kb, vmsize_kb);

    
    /* Firstly load hashes to the dictionary */
    
    
    dict* d = dict_new(NSLOTS_MY_NODE); /* This is done by python script */

    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    printf("recv #%d memory usage after dict creation: ram:"
	   "%ld kb vm: %ld kb\n",
	   myrank, vmrss_kb, vmsize_kb);

    double start_time = wtime();
    char file_name[FILE_NAME_MAX_LENGTH];
    snprintf(file_name, sizeof(file_name), "data/receive/digests/%d", myrank);
    FILE* fp = fopen(file_name, "r");
    load_file_to_dict(d, fp);


    printf("recv #%d dict read in %0.2fsec. It has %lu elms, file has %lu elms."
	   " dict.nslots = %lu, filling rate=%f \n",
	   myrank, wtime() - start_time, d->nelements,
	   d->nelements_asked_to_be_inserted,
	   d->nslots, ((float) d->nelements)/d->nslots);
    
    fclose(fp);

    
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    printf("reciver #%d memory usage after load: ram: %ld kb vm: %ld kb\n"
	   "============================================================\n",
	   myrank, vmrss_kb, vmsize_kb);
    // 000000000000000000000000000000000000000000000000000000000000000000000000000
    //-------------------------------------------------------------------------+
    // I'm a receiving process: receive hashes, probe them, and send candidates 
    // Process Numbers: [0,  NSERVERS - 1]
    //-------------------------------------------------------------------------+
    // Receive templates from all senders only once:
    
    u8* templates = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nproc_snd);
    receiver_process_get_template(myrank, nproc, nproc_snd, templates);


    /* listen to sender and save candidates */
    receiver_process_task(d, myrank, nproc, nproc_snd, templates);
    

    /* send the messages to the archive! */
    /* printf("I am receiver with rank %d is going to send to archive\n", myrank); */
    /* send_candidates_to_archive(myrank); */
    /* printf("I am receiver rank %d is done sending to archive\n", myrank); */


    free(templates);
  }



  // The end that will never be reached in any case!
  printf("process #%d reached the end (should be a receiver)\n", myrank);
  MPI_Finalize();
  return 0;
}


// ------------------- Auxililary functions phase iii --------------------------

/* /\* 1 if  dgst1 > dgst2, -1 if dgst1<dgist2, 0 if dgst1==dgst2 *\/ */
/* int cmp_dgst(void const* dgst1, void const* dgst2){ */
/*   return memcmp(dgst1, dgst2, N); /\* comparison order: low bytes first *\/ */
/* } */

/* /\* return index of key if it is found, -1 otherwise*\/ */
/* int64_t linear_search(u8 *key, u8 *array, size_t array_len, size_t key_len) */
/* { */
/*   for (size_t i=0; i<array_len; i += key_len) { */
/*     if ( 0 == memcmp(key, &array[i], key_len) ) */
/*       return i; */
/*   } */
  
/*   return -1; /\* not found *\/ */
/* } */


/* void print_byte_array(u8* array, size_t nbytes) */
/* { */
/*   for (size_t i=0; i<nbytes; ++i)  */
/*     printf("0x%02x, ",  array[i]); */
/*   puts(""); */
/* } */





