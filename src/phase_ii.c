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



//---------------------------- UTILITY FUNCTIONS -------------------------------


static inline int lookup_multi_save(dict *d,
				    u8 *stream,
				    u8 init_message[HASH_INPUT_SIZE],
				    size_t npairs,
				    FILE *fp)
{ // Given stream={msg1||dgst1,..} from specific rank, probes each dgst, if   |
  // prope(dgst)=/=0 store its related msg in fp, and returns the number of   |
  // stored messages.                                                         |
  // -------------------------------------------------------------------------+
  // INPUT:                                                                   |
  // - dict                                                                   |
  // - stream                                                                 |
  // - npairs: how many pairs (msg, dgst) in the stream.                      |
  // -------------------------------------------------------------------------+
  // NOTES:                                                                   |
  // |msg| = NWORDS_INPUT*WORD_SIZE (config.h)                                |
  // |dgst| = N-DEFINED_BYTES (should be N,but we skip know bits .e.g nserver)|
  // -------------------------------------------------------------------------+
  static int one_pair_size = HASH_INPUT_SIZE + (N-DEFINED_BYTES);
  static int msg_size = HASH_INPUT_SIZE;
  
  /* how many messages give positive ans when their dgst gets probes */
  int npositive_msgs = 0;
  int tmp = 0;
  for (size_t i=0; i<npairs; ++i){
    /* dictionary only read |dgst| bytes by default  */
    tmp =  dict_has_elm(d, stream+i*one_pair_size+msg_size);

    if (tmp){ /* positive probe */
      // reconstruct the message:
      /* I feel we should pass this as an argument */
      static char M[HASH_INPUT_SIZE];
      memcpy(M, init_message, HASH_INPUT_SIZE);
      /* set the counter part */
      memcpy(M, stream+i*one_pair_size, sizeof(CTR_TYPE));
      /* finally write the reconstructed message */
      fwrite(M, sizeof(u8), msg_size, fp);
      fflush(fp); /* ensures it's written */
      
      ++npositive_msgs;
    }
  }
  return npositive_msgs;
}

static void find_hash_distinguished(u8 M[HASH_INPUT_SIZE], /* in, out*/
				    WORD_TYPE Mstate[NWORDS_STATE], /* out*/
				    const size_t dist_test /* in */)
{
  // ==========================================================================+
  // Summary: Generates hashes till a distinguished point is found. Mutates M  |
  //          and resests Mstate in the process. Suitable for phase ii use only|
  // --------------------------------------------------------------------------+
  // We start with message M, computed hash(M) if it is not a distinguished    |
  // point then we change M slightly. Increment the first 64 bits of M by 1    |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // - M[16] : initial value of the random message, we hash this in the first  |
  //           in first trial, then change it to M xor ctr                     |
  // - Mstate: This should be the initial state of sha256.                     |
  // - dist_test: (2^nzeros - 1), e.g. a point is distinguished if it has 3    |
  //              at the end, then the test is 2^3 - 1 = 7                     |
  // --------------------------------------------------------------------------+
  // WARNING: this function can't deal with more than 32zeros as dist_test     |
  // NOTE: we are only working on the first word                               |
  // --------------------------------------------------------------------------+
  // TODO: use hash_multiple instead                                           |
  // --------------------------------------------------------------------------+

  /* no need to construct init state with each call of the function */ 
  const static WORD_TYPE init_state[NWORDS_STATE] = {HASH_INIT_STATE};
  


  /* increments the first sizeof(CTR_TYPE)*8 bits of M by 1 */
  CTR_TYPE* ctr_pt = (CTR_TYPE*) M;
  
  while (1) { /* loop till a dist pt found */
    ++(ctr_pt[0]);
   
    memcpy(Mstate, init_state, 32);
    /* todo  use hash multiple */
    /* figure out the number of words from config.h */
    hash_single(Mstate, (u8*) M);
    
    /* is its digest is a distinguished pt? */
    /* see 1st assumption in config.h  */
    if ( (((WORD_TYPE*) Mstate)[0] & dist_test) == 0){

      return; /* we got our distinguished digest */
    }

      
  }
}

// -----------------------------------------------------------------------------

static inline u32 to_which_server(u8 MState[NWORDS_DIGEST*WORD_SIZE])
{
  // ==========================================================================+
  // Summary: Given a state. Decide to which server we send it to.             |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `Mstate` : Array of bytes of the digest.                                  |
  // --------------------------------------------------------------------------+
  // Given a state. Decide to which server we send to                          |
  // """ One potential idea is to do:                                          |
  // server <-- h % nserver                                                    |
  // h'     <-- h / nserver """                                                |
  //
  // --------------------------------------------------------------------------+
  
  const static u32 ones_nservers = (1LL<<LOG2_NSERVERS) - 1;
  /* 1- convert MState to a WORD_TYPE (dist_bits || nserver || rest)32bits */
  /* 2- remove the distinguished bits by shifting  (nserver || rest ) */
  /* 3- keep only the bits that holds nserver (nserver) it may have extra bit */
  /* 4- Compute server number by taking computing mod nservers */
  u32 snd_to_server  = ( (((WORD_TYPE*)MState)[0] >> DIFFICULTY)
			 & ones_nservers) % NSERVERS;

  return snd_to_server;
}
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

  // @todo address sanitizer detected stack buffer overflow, stream_pt[6], L_IN_BYTES=3, VAL_SIZE_BYTES=4
  u8 stream_pt[N-DEFINED_BYTES]; /* @todo I believe the issue is in here */
  /* add as many hashes as possible */ 
  while ( !feof(fp) ){
    fread(stream_pt, sizeof(u8), N-DEFINED_BYTES, fp);
    /* it adds the hash iff nprobes <= NPROBES_MAX */
    dict_add_element_to(d, stream_pt); 
  }
  
  // fclose(fp); // don't close the file 
}
// -----------------------------------------------------------------

void send_random_message_template(u8 M[HASH_INPUT_SIZE])
{ /* Send M immediately and clear the memory at the end  */

  MPI_Request requests[NSERVERS]; /* they will be used only here  */
  MPI_Status statuses[NSERVERS];
  // how about collective communications? I can't find a simple way using them.
  for (int i=0; i<NSERVERS; ++i) {
    MPI_Isend(M, HASH_INPUT_SIZE, MPI_UNSIGNED_CHAR, i,
	      TAG_RANDOM_MESSAGE, MPI_COMM_WORLD, &requests[i]);
  }
  MPI_Waitall(NSERVERS, requests, statuses);
} /* clear stack variables */



void sender_process_task(u8 M[HASH_INPUT_SIZE], int myrank)
{
  /* generate hashes and send them to a server as soon its buffer is complete */
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}
  // dgst := N bytes of the digest
  // ctr  := (this is a shortcut that allows us not to send the whole message)
  //---------------------------------------------------------------------------+
  WORD_TYPE Mstate[NWORDS_STATE] = {HASH_INIT_STATE};
  int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                    + sizeof(CTR_TYPE); /* |dgst| + |ctr| - |known bits|*/

  int server_number;


  // ---- Part1 : initializing sending buffers 
  /* int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES) */
  /*                      + sizeof(CTR_TYPE); /\* |dgst| + |ctr| *\/ */

  u8* snd_buf = (u8*) malloc(one_pair_size
			     *PROCESS_QUOTA
			     *NSERVERS);

  /* Init attaching buffers for the buffered send */
  int buf_attached_size = one_pair_size*PROCESS_QUOTA + MPI_BSEND_OVERHEAD;
  u8* buf_attached = (u8*) malloc(buf_attached_size);
  MPI_Buffer_attach(buf_attached, buf_attached_size);   

  /* Decide where to place the kth digest in server i buffer */
  size_t offset = 0;
  const u64 mask_test = (1LL<<DIFFICULTY) - 1;

  /* pos i: How many messages we've generated to be sent to server i? */
  u32* servers_ctr = (u32*) malloc(sizeof(u32)*NSERVERS);
  /* set number of generated messages before sending to 0 */
  memset(servers_ctr, 0, sizeof(u32)*NSERVERS); 

  //-------------------------------------------------------------------------+
  //------- Part 2 : Generate hashes and send them  -------
  //-------------------------------------------------------------------------+
 
  printf("send #%d: done init mpi, now going to generate hashes \n", myrank);
  
  // find_hash_distinguished_init(); 
  int snd_ctr = 0;


  while(1) { /* when do we break? never! */
    /* Find a message that produces distinguished point */

    find_hash_distinguished( M, Mstate, mask_test);

    //+ decide to which server to add to? 
    server_number = to_which_server((u8*) Mstate);

    /* 1st term: go to server booked memory, 2nd: location of 1st free place*/
    offset = server_number*one_pair_size + servers_ctr[server_number];
    // recall que one_pair_size =  |dgst| + |ctr| - |known bits|
    
    /* save N-DEFINED_BYTES of MState in: snd_buf_dgst[offset] */
    memcpy( ((u8*)Mstate) + DEFINED_BYTES, /* skip defined bytes */
	     &snd_buf[offset], /* copy digest to snd_buf[offset] */
	     N-DEFINED_BYTES );

    /* record the counter, above we've recorded N-DEFINED_BYTES  */
    memcpy(M,
	   &snd_buf[ offset + (N-DEFINED_BYTES) ],
	   sizeof(CTR_TYPE) );

    /* this server has one more digest */
    ++servers_ctr[server_number];

    
    if (servers_ctr[server_number] == PROCESS_QUOTA){
      ++snd_ctr;
      
      /* printf("rank %d: sending to %d, snd_ctr=%d, it took %fsec\n", */
      /* 	     myrank, server_number, snd_ctr, wtime() - time_start); */
      /* time_start = wtime(); */
      /* we have enough messages to send to server (server_number) */
      MPI_Send(snd_buf,
		PROCESS_QUOTA*one_pair_size,
		MPI_UNSIGNED_CHAR,
		server_number,
		TAG_SND_DGST,
		MPI_COMM_WORLD);
      /* printf("rank %d: sending done to %d, snd_ctr=%d, it took %fsec\n", */
      /* 	     myrank, server_number, snd_ctr, wtime() - time_start); */
      /* time_start = wtime(); */
      // todo do we need here buffer detach? 
      /* It is enough to reset the counter. The memroy will be rewritten */  
      servers_ctr[server_number] = 0;
    }
  } /* todo why do I stop? it is not specified!pp */

  MPI_Buffer_detach(buf_attached, &buf_attached_size);

  free(snd_buf);
  free(buf_attached);
  free(servers_ctr);
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
    printf("receiver %d is listening to %d for the template\n", myrank, i + NSERVERS);
    MPI_Recv(&templates[i*HASH_INPUT_SIZE],
	     HASH_INPUT_SIZE,
	     MPI_UNSIGNED_CHAR,
	     i + NSERVERS,
	     TAG_RANDOM_MESSAGE,
	     MPI_COMM_WORLD,
	     &status);

    
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
  
  char file_name[128]; /* "data/send/messages/%d" */
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
  int sender_name = 0;


  
  
  long vmrss_kb, vmsize_kb;
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  printf("reciver #%d memory usage after mpi init: ram: %ld kb vm: %ld kb\n",
	 myrank, vmrss_kb, vmsize_kb);
  




  //---------------------------------------------------------------------------+
  // --- Part 3: Receive digests and probe them
  //---------------------------------------------------------------------------+
  printf("recv#%d is listening\n", myrank);
  // listen the first time
  MPI_Recv(rcv_buf,
	   rcv_array_size,
	   MPI_UNSIGNED_CHAR,
	   MPI_ANY_SOURCE,
	   TAG_SND_DGST,
	   MPI_COMM_WORLD,
	   &status);
  // copy the received message and listen immediately
  memcpy(lookup_buf, rcv_buf, rcv_array_size);
  sender_name = status.MPI_SOURCE; // who sent the message?


  while (NNEEDED_CND_THIS_SERVER  > nfound_cnd) {

    //+ receive messages from different processors
    MPI_Irecv(rcv_buf, /* store in this location */
	      rcv_array_size, 
	      MPI_UNSIGNED_CHAR,
	      MPI_ANY_SOURCE, /* sender */
	      TAG_SND_DGST, 
	      MPI_COMM_WORLD,
	      &request);

    //+ probe these messages and update the founded candidates
    printf("recv#%d is going to probe the dict\n", myrank);
    nfound_cnd += lookup_multi_save(d, /* dictionary to look inside */
				    rcv_buf, /* messages to search in d */
				    &templates[sender_name],
				    PROCESS_QUOTA,/* how many msgs in rcv_buf */
				    fp); /* file to record cadidates */

    if (nfound_cnd - old_nfound_candidates > 0) {
      printf("receiver #%d found %lu candidates\n", myrank, nfound_cnd);
      old_nfound_candidates = nfound_cnd;
    }
    MPI_Wait(&request, &status);
    sender_name = status.MPI_SOURCE; // get the name of the new sender
    memcpy(lookup_buf, rcv_buf, rcv_array_size);
    }


  // good job
  free(rcv_buf);
  free(indices);
  fclose(fp);
  printf("recv #%d done a good job\n", myrank);
}






// 
// -----------------------------------------------------------------------------


void phase_ii()
{ 
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
  
  
  /* send digests, 1 digest ≡ 256 bit */
  // int nthreads = omp_get_max_threads(); // no need for this?!

  // Who am I? sender, receiver, or archive.
  if (myrank >= NSERVERS){
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
    printf("sender %d got template\n", myrank);
    print_char(M,HASH_INPUT_SIZE);

    
    /* Send the initial input to all receiving servers */
    printf("sender #%d is going to send its template\n", myrank);
    send_random_message_template(M); 
    printf("sender #%d done sending template!\n", myrank);
    /* generate hashes and send them to a servers */

    sender_process_task(M, myrank); /* never ends :) */
  }

  //+ todo receive processors
  //+ load dgsts from file to dictionary

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
    char file_name[40];
    snprintf(file_name, 40, "data/receive/digests/%d", myrank);
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
    //-------------------------------------------------------------------------+
    // I'm a receiving process: receive hashes, probe them, and send candidates 
    // Process Numbers: [0,  NSERVERS - 1]
    //-------------------------------------------------------------------------+
    // Receive templates from all senders only once:
    
    u8* templates = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nproc_snd);
    receiver_process_get_template(myrank, nproc, nproc_snd, templates);
    

    /* listen to sender and save candidates */
    printf("I am reciever with rank %d i am listening \n", myrank);
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
  return;
}


// ------------------- Auxililary functions phase iii --------------------------

/* 1 if  dgst1 > dgst2, -1 if dgst1<dgist2, 0 if dgst1==dgst2 */
int cmp_dgst(void const* dgst1, void const* dgst2){
  return memcmp(dgst1, dgst2, N); /* comparison order: low bytes first */
}

/* return index of key if it is found, -1 otherwise*/
size_t linear_search(u8 *key, u8 *array, size_t array_len, size_t key_len)
{
  for (size_t i=0; i<array_len; i += key_len) {
    if ( 0 == memcmp(key, &array[i], key_len) )
      return i;
  }
  
  return -1; /* not found */
}


void print_byte_array(u8* array, size_t nbytes)
{
  for (size_t i=0; i<nbytes; ++i) 
    printf("0x%02x, ",  array[i]);
  puts("");
}




int main() {  phase_ii(); }
// @todo get the number of stored candidates at each run

