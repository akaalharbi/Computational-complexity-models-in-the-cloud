// phase ii: high level overview 
// four types of processors: senders (the majority), receivers (#NSERVERS),
// and archive (1).
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
// `archive`: -(listen) stay posted for a message to be received from `receiver`
//             when receiving a message, write them immediately to a file.





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
    /* if ( ((WORD_TYPE*) Mstate)[0] && dist_test == 0){ */
    /*   puts("found a dist point"); */
    /*   return; /\* we got our distinguished digest *\/ */
    /* } */
    return;
      
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

  u8 stream_pt[N-DEFINED_BYTES];
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
  int buf_attached_size = NSERVERS
                        *(one_pair_size*PROCESS_QUOTA + MPI_BSEND_OVERHEAD);
  
  u8* buf_attached = (u8*) malloc(buf_attached_size);
    

  /* Decide where to place the kth digest in server i buffer */
  size_t offset = 0;
  // h is distinguished iff (h & mask_test) == 0 
  const u64 mask_test = (1LL<<DIFFICULTY) - 1;

  /* pos i: How many messages we've generated to be sent to server i? */
  u32* servers_ctr = (u32*) malloc(sizeof(u32)*NSERVERS);
  /* set number of generated messages before sending to 0 */
  memset(servers_ctr, 0, sizeof(u32)*NSERVERS); 

  //-------------------------------------------------------------------------+
  //------- Part 2 : Generate hashes and send them  -------
  //-------------------------------------------------------------------------+
  MPI_Buffer_attach(buf_attached, buf_attached_size);
  printf("rank %d: done init mpi, now going to generate hashes \n", myrank);
  // find_hash_distinguished_init(); 
  while(1) { /* when do we break? never! */
    /* Find a message that produces distinguished point */
    find_hash_distinguished( M, Mstate, mask_test);

    //+ decide to which server to add to? 
    server_number = to_which_server((u8*) Mstate);

    /* 1st term: go to server booked memory, 2nd: location of 1st free place*/
    offset = server_number*one_pair_size + servers_ctr[server_number];
    // recall que one_pair_size =  |dgst| + |ctr| - |known bits|
    
    /* save N-DEFINED_BYTES of MState in: snd_buf_dgst[offset] */
    memcpy( ( (u8*)Mstate) + DEFINED_BYTES, /* skip defined bytes */
	      (snd_buf+offset), /* copy digest to snd_buf[offset] */
	      N-DEFINED_BYTES );

    /* record the counter, above we've recorded N-DEFINED_BYTES  */
    memcpy(M,
	   ( (snd_buf+offset) + (N-DEFINED_BYTES) ),
	   sizeof(CTR_TYPE) );

    /* this server has one more digest */
    ++servers_ctr[server_number];


    
    if (servers_ctr[server_number] == PROCESS_QUOTA){
      printf("rank %d: sending to %d\n", myrank, server_number);
      /* we have enough messages to send to server (server_number) */
      MPI_Bsend(snd_buf,
		PROCESS_QUOTA*one_pair_size,
		MPI_UNSIGNED_CHAR,
		server_number,
		TAG_SND_DGST,
		MPI_COMM_WORLD);
               
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




void receiver_process_task(dict* d, int myrank, int nproc, int nproc_snd )
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
  
  char file_name[25]; /* "data/send/messages/%d" */
  snprintf(file_name, sizeof(file_name), "data/send/messages/%d", myrank );
  FILE* fp = fopen(file_name, "w");
  int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                    + sizeof(CTR_TYPE); /* |dgst| + |ctr| */
  
  //+ add initial random message buffers
  u8* rcv_buf = (u8*) malloc(one_pair_size
			     *nproc_snd /* we've more senders */
			     *PROCESS_QUOTA);
    

  MPI_Status* statuses = (MPI_Status*)malloc(sizeof(MPI_Status)*nproc_snd);
  MPI_Request* requests = (MPI_Request*)malloc(sizeof(MPI_Request)*nproc_snd);
  int* indices = (int*)malloc(sizeof(int)*nproc);
  size_t idx = 0;
  size_t buf_idx = 0;
  size_t msg_idx = 0;
  int flag = 0; /* decide if all messages have been received */
  int outcount = 0; /* MPI_Waitsome  how many buffers have we received so far */
  size_t nfound_cnd = 0;
  size_t rcv_array_size = one_pair_size*PROCESS_QUOTA;
  u8* initial_inputs = (u8*) malloc(sizeof(u8)*HASH_INPUT_SIZE*nproc_snd);

  long vmrss_kb, vmsize_kb;
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  printf("reciver #%d memory usage after mpi init: ram: %ld kb vm: %ld kb\n",myrank, vmrss_kb, vmsize_kb);

  printf("receiver %d done with initialization for mpi\n", myrank);
  
  //---------------------------------------------------------------------------+
  // --- part 2: receive the initial inputs from all generating processors 
  //---------------------------------------------------------------------------+

  // process : 0 -> NSERVERS - 1 (receivers),
  //         : NSERVERS (archive)
  //         : NSERVERS + 1 -> nproc (senders)
  for (int i=0; i<nproc_snd; ++i){
    printf("receiver %d is listening to %d\n", myrank, i + NSERVERS + 1);
    MPI_Recv(&initial_inputs[i*HASH_INPUT_SIZE],
	     HASH_INPUT_SIZE,
	     MPI_UNSIGNED_CHAR,
	     i + NSERVERS + 1,
	     TAG_RANDOM_MESSAGE,
	     MPI_COMM_WORLD,
	     statuses);
    printf("recv #%d has template of sender #%d \n", myrank, i + NSERVERS + 1);
    print_char(&initial_inputs[i*HASH_INPUT_SIZE], HASH_INPUT_SIZE);
  }
  printf("receiver %d received all templates\n", myrank);
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  printf("reciver #%d memory usage after receiving temp: ram: %ld kb vm: %ld kb\n",myrank, vmrss_kb, vmsize_kb);

  //---------------------------------------------------------------------------+
  // --- Part 3: Receive digests and probe them
  //---------------------------------------------------------------------------+
  // enclose what's below within while loop to receive multiple times
  // Warning: not correct! nfound_collisions should be shared between all
  while (NNEEDED_CND_THIS_SERVER  > nfound_cnd) {
    //+ receive messages from different processors
  
    for (int i = 0; i<nproc_snd; ++i) {
      MPI_Irecv(&rcv_buf[i*rcv_array_size], /* store in this location */
		rcv_array_size, 
		MPI_UNSIGNED_CHAR,
		i + NSERVERS + 1, /* sender */
		TAG_SND_DGST, 
		MPI_COMM_WORLD,
	       &requests[i]);
    }
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    // printf("reciver #%d memory usage after irecv: ram: %ld kb vm: %ld kb\n",myrank, vmrss_kb, vmsize_kb);

    //+ probe these messages 
    while (!flag) {/*receive from any server and immediately treat their dgst*/
      MPI_Waitsome(nproc_snd, requests, &outcount, indices, statuses);
      /* check if ther is no further messages */
      MPI_Testall(nproc_snd, requests, &flag, statuses);
      
      // printf("receiver %d I received from %d senders", myrank, outcount);

      for (int j = 0; j<outcount; ++j){
	printf("receiver %d going to treat %d", myrank, j);
	idx = indices[j]; /* number of sending process */
	buf_idx = idx*one_pair_size*PROCESS_QUOTA;
	msg_idx = idx*HASH_INPUT_SIZE; // location of message template of sender idx
	/* treat received messages that are completed */
	//+ probe dictionary
	//+ if it is positive:
	//+ reconstruct the messages
	//+ record message and digest in file fp
	nfound_cnd += lookup_multi_save(d,
					&rcv_buf[buf_idx],
					&initial_inputs[msg_idx],
					PROCESS_QUOTA,
					fp);
      } // end treating received messages 
    } // flag == 1 when buffers have been received, thus call mpi_irecv again
  } // end found enough collisions

  // good job

  free(rcv_buf);
  free(statuses);
  free(requests);
  free(initial_inputs);
  free(indices);
  fclose(fp);
}



void send_candidates_to_archive(int myrank)
{
  /* allocate buffer for messages that found s.t. their dgst give postivie ans*/
  
  size_t snd_buf_size = sizeof(u8)*HASH_INPUT_SIZE*NNEEDED_CND_THIS_SERVER;
  u8* snd_buf = (u8*) malloc(snd_buf_size);

  char msg_file_name[40];
  /* file name := myrank_t(time_stamp) */
  snprintf(msg_file_name, sizeof(msg_file_name),
	   "data/send/messages/%d_t%llu", myrank, (u64) wtime());
  
  FILE* fp = fopen(msg_file_name, "r");
  fread(snd_buf, snd_buf_size, 1, fp); /* should matches with the file size */
  MPI_Send(snd_buf, snd_buf_size, MPI_UNSIGNED_CHAR,
	   ARCHIVE_SERVER,
	   TAG_MESSAGES_CANDIDATES,
	   MPI_COMM_WORLD);
  fclose(fp);
  
  return;

}


void archive_receive()
{

    char archive_file_name[] = "data/receive/messages/archive";
    int actual_rcv_count = 0;
    
    /* MPI_Status status; */
    MPI_Status* statuses = (MPI_Status*)malloc(sizeof(MPI_Status)*NSERVERS);
    MPI_Request* requests = (MPI_Request*)malloc(sizeof(MPI_Request)*NSERVERS);
    int* indices = (int*)malloc(sizeof(int)*NSERVERS);

    int flag = 0; /* decide if all messages have been received */
    int outcount = 0; /* MPI_Waitsome  how many buffers have we received so far */

    /* maximum possible receiving message  */
    u8* rcv_buf = malloc(sizeof(u8)*MAX_CND_PER_SERVER*NSERVERS);
    FILE* fp = fopen(archive_file_name, "a");
    
    for (int i = 0; i<NSERVERS; ++i) {
      MPI_Irecv(rcv_buf+i*MAX_CND_PER_SERVER,
		MAX_CND_PER_SERVER,
		MPI_UNSIGNED_CHAR,
		i,
		TAG_MESSAGES_CANDIDATES,
		MPI_COMM_WORLD,
		&requests[i]);
    }

    /* MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &rcv_count); */
    /* fwrite(rcv_buf, rcv_count, 1, fp); */
    /* fclose(fp); */
    
    while (!flag) {/*receive from any server and immediately treat their dgst*/
      MPI_Waitsome(NSERVERS, requests, &outcount, indices, statuses);
      MPI_Testall(NSERVERS, requests, &flag, statuses);
      for (int j = 0; j<outcount; ++j){
	/* How many messages have we received? */
	MPI_Get_count(&statuses[indices[j]],
		      MPI_UNSIGNED_CHAR,
		      &actual_rcv_count);

	/* Number of received messages */
	actual_rcv_count = actual_rcv_count / HASH_INPUT_SIZE;
	fwrite(&rcv_buf[indices[j]],
	       actual_rcv_count*HASH_INPUT_SIZE,
	       1,
	       fp);

        fflush(fp);
      } // end treating received messages 
    } // flag == 1 when buffers have been received, thus call mpi_irecv again
    fclose(fp);
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
  int nproc_snd = nproc - NSERVERS - 1;
  
  
  /* send digests, 1 digest ≡ 256 bit */
  // int nthreads = omp_get_max_threads(); // no need for this?!

  // Who am I? sender, receiver, or archive.
  if (myrank > NSERVERS){
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
    printf("I am sender with rank %d i am going to send my template \n", myrank);
    send_random_message_template(M); 

    /* generate hashes and send them to a servers */
    printf("I am sender with rank %d i am going to generate hashes  \n", myrank);
    sender_process_task(M, myrank); /* never ends :) */
  }

  //+ todo receive processors
  //+ load dgsts from file to dictionary

  while (myrank < NSERVERS){ /* receiver, repeat infinitely  */
    long vmrss_kb, vmsize_kb;
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);

    /* Firstly load hashes to the dictionary */
    printf("reciver #%d memory usage beginning: ram: %ld kb vm: %ld kb\n",myrank, vmrss_kb, vmsize_kb);
    dict* d = dict_new(NSLOTS_MY_NODE); /* This is done by python script */
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    printf("reciver #%d memory usage after dict creation: ram: %ld kb vm: %ld kb\n", myrank, vmrss_kb, vmsize_kb);
    char file_name[40];
    snprintf(file_name, 40, "data/receive/digests/%d", myrank);
    FILE* fp = fopen(file_name, "r");
    load_file_to_dict(d, fp);
    fclose(fp);

    
    get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
    printf("reciver #%d memory usage after load: ram: %ld kb vm: %ld kb\n", myrank, vmrss_kb, vmsize_kb);
    //-------------------------------------------------------------------------+
    // I'm a receiving process: receive hashes, probe them, and send candidates 
    // Process Numbers: [0,  NSERVERS - 1]
    //-------------------------------------------------------------------------+

    
    /* Listen to senders till we accumulated enough candidates */
    printf("I am reciever with rank %d i am listening \n", myrank);
    receiver_process_task(d, myrank, nproc, nproc_snd);
    
    printf("I am receiver with rank %d is going to send to archive\n", myrank);

    /* send the messages to the archive! */
    send_candidates_to_archive(myrank);
    printf("I am receiver rank %d is done sending to archive\n", myrank);
  }

  while (myrank == NSERVERS){ /* I'm THE archive*/
    //+ todo truncate the messages file if it not multiple to HASH_INPUT_SIZE
    //+ this happen if the server was shut down during writing
    printf("I am the archive listening \n");
    archive_receive();
    printf("I am archive registered messages on disk\n");
  }


  // The end that will never be reached in any case!
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
