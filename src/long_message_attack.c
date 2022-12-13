// Long message attack





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


/* // MPI sending and receiving tags */
/* #define TAG_DICT_SND 0 /\* send digest tag *\/ */
/* #define TAG_RANDOM_MESSAGE 1 */
/* #define TAG_SND_DGST 2 */
/* #define TAG_MESSAGES_CANDIDATES 3 */

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
  // TODO: Don't use hash_multiple instead                                     |
  // --------------------------------------------------------------------------+

  /* no need to construct init state with each call of the function */ 
  const static WORD_TYPE init_state[NWORDS_STATE] = {HASH_INIT_STATE};
  


  /* increments the first sizeof(CTR_TYPE)*8 bits of M by 1 */
  CTR_TYPE* ctr_pt = (CTR_TYPE*) M;
  
  while (1) { /* loop till a dist pt found */
    ++(ctr_pt[0]);
   
    memcpy(Mstate, init_state, 32);
    /* todo  */
    hash_single(Mstate, (u8*) M);
    
    /* is its digest is a distinguished pt? */
    /* see 1st assumption in config.h  */
    if ( ((u32*) Mstate)[0] && dist_test == 0)
      return; /* we got our distinguished digest */
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



void phase_i_store(size_t server_capacity[]){
		   /* const size_t n, */
		   /* size_t global_difficulty, */
		   /* size_t nservers */
                   /* parameters above moved to config.h */
  // ==========================================================================+
  // Summary:                                                                  |
  // Hash a long message of zeros. Store the digest, h, in a file k where      |
  // 0<= k < nservers. It will store N bytes of digest (N defined in config.h) |
  // To decide which server gets the digest h, compute k := h1 mod  nservers   |
  // where h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest.     |
  // --------------------------------------------------------------------------+
  // INPUTS: CAPITAL LETTERS input are defined in config.h                     |
  //                                                                           |
  // `server_capacity[]` : array of size nservers,  entry i contains how many  |
  //                       blocks server i will store in its dictionary        |
  // `DIFFICULTY` : Number of bits that are 0 in the first word A:=state[0]    |
  //                       i.e. is state[0]&(2**global_difficulty - 1) == 0?   |
  // `NSERVERS` : how many servers we should prepare for                       |
  // --------------------------------------------------------------------------+
  // NOTE: Bits corresponding to distinguished point and server number are not |
  //       stored in the file.                                                 |
  // NOTE: endinaness: u32 A[2] = {x, y} then (uint64*) A = { (y<<32) | x) }   |
  // --------------------------------------------------------------------------+
  // TODO:                                                                     |
  // ==========================================================================+
  



  /// ----------------------- INIT ---------------------------///
  /// 1- INIT numerical and bytes variables:
  /* phase will rehash the long message again in parallel */
  /* here we define how many parllel processors in phase iii */
  size_t ncores = 14; 
  size_t k =  0; // server index
  int should_NOT_stop = 1;
  size_t nhashes_stored = 0; // 

  u32 ones = (1LL<<DIFFICULTY) - 1;



  /* Actually we have discussed that we can fix the number of nhashes per */
  /* Eventhough they have different capacities, since we allow server to  */
  /* discard excessive hasing */
  size_t nhashes=0;
  for (size_t i=0; i<NSERVERS; ++i) {
    /* this should coincide with number of hashes in config.h */
    nhashes += server_capacity[i]; 
  }

  /* record the whole state after each each interval has passed */
  size_t interval = nhashes>>10; 
  
  /// timing variables
  double start = 0;
  double elapsed = 0;

  // INIT SHA256 
  u8 M[HASH_INPUT_SIZE] = {0};

  // store the hash value in this variable
  WORD_TYPE state[NWORDS_STATE] = {HASH_INIT_STATE};

  /* treat state as bytes (have you considered union?) */
  u8* stream_pt = (u8*) state; 

  /// INIT FILES: server files that will be send later and state
  char file_name[40]; // more than enough to store file name
  char states_file_name[40];
  
  FILE* data_to_servers[NSERVERS];
  FILE* states_file;
 
  // TOUCH FILES ON THE DISK
 
  fopen("data/states", "w");
  /* create a string that will become a file name  */
  snprintf(states_file_name,
	   sizeof(states_file_name),
	   "data/%llu_state",
	   (u64) wtime());

  states_file = fopen(states_file_name, "w");
  /* fclose(states_file); // we will open this file again in few occasions */
  
  for (size_t i=0; i<NSERVERS; ++i) {
    //edit file name according to server i
    snprintf(file_name, sizeof(file_name), "data/upload/%lu", i);
    printf("file_name=%s\n", file_name);

    data_to_servers[i] = fopen(file_name, "w");
    nhashes_stored += server_capacity[i];
  }

  // Init coutners before the beginning of the attack
  interval = nhashes_stored / ncores;
  printf("interval=%ld, nhashes_stores=%ld, ncores=%ld\n",
	 interval, nhashes_stored, ncores);
  nhashes_stored = 0; // we have not recorded any hash yet




  /// ----------------- PHASE I: Part 1/2   ------------------------  ///
  // First phase hash an extremely long message
  // M0 M1 ... M_{2^l}, Mi entry will evaluated on the fly
  // Store the hashes in file correspond to some server k
  start = wtime();
  /* if one server gets filled, it will */
  while (should_NOT_stop) {
    // hash and extract n bits of the digest
    hash_single(state, M);

    if ((state[0] & ones) == 0){ /* it is a distinguished point */
      
      /* Decide which server is responsible for storing this digest */
      k = to_which_server((u8*) state);
	//( (state[0]>>DIFFICULTY) & ones_nservers) % NSERVERS;

      /* Recall that: */
      /* h = (dist_pt) || h1:=b0 ... b_ceil(log2(nservers)) || the rest   */
      fwrite(stream_pt+DEFINED_BYTES, /* start  from "the rest" see above */
	     sizeof(u8), /* smallest moving unit */
	     N-DEFINED_BYTES, /* len( (dist_pt)|| h1 ) = DEFINED_BITS */
	     data_to_servers[k]);
      
      server_capacity[k] -= 1; // We need to store less blocks now
      ++nhashes_stored;

      
      // decide should not stop or should stop?
      /* not the most optimal implementation */
      should_NOT_stop = 0;
      for (size_t i=0; i<NSERVERS; ++i) {
	/*--------------------------------------------------------------------*/
	/* Since we reduce the server capacity by one each time when we add   */
	/* to its file. If any server has capacity larger than zero then it   */
	/* means that we should. continue hashing till all servers capacities */
	/* have been rached.                                                  */
	/*------------------------------------------------------------------- */
	/* should_NOT_stop == 0 iff all servers_capacities are 0; */
	should_NOT_stop |= (server_capacity[i] > 0) ;
      }
    

      // + save states after required amount of intervals
      
      if (nhashes_stored % interval == 0) {
	/* FILE* states_file = fopen(states_file_name, "a"); */
	
	/* Record the whole state */
	fwrite((WORD_TYPE*) stream_pt, sizeof(WORD_TYPE), NWORDS_STATE, states_file);


	// We would like to flush the data disk as soon we have them
	fflush(states_file);
      }
    }
    
  }
  fclose(states_file);


  elapsed = wtime() - start;
  /// write it in a file  ///
  printf("done in %fsec, ", elapsed);
}


void phase_i_load(dict *d, size_t fp_nhashes, FILE *fp)
{
  // ==========================================================================+
  // Summary: Load hashes from the file *fp, and store them in dict *d         |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d`: Dictionary that will keep elements from *fp.                        |
  // `*fp` : File contain number of hashes larger than nelements               |
  // `fp_nhashes`  : The nhashes available in the file.                        |
  // `mem_nhashes` : The max number of hashes we are allowed to store in this  |
  //                   machine. The dictionary *d may not accepts all elements.|
  // ==========================================================================+


  /* Check that file exists, the file comes from external resources */  
  if (!fp){
    puts("I have been given a file that lives in nowhere");
    return; // raise an error instead
  }

  u8 stream_pt[N-DEFINED_BYTES];
  /* add as many hashes as possible */
  for (;; fp_nhashes--){
    fread(stream_pt, sizeof(u8), N-DEFINED_BYTES, fp);
    /* it adds the hash iff nprobes <= NPROBES_MAX */
    dict_add_element_to(d, stream_pt); 
  }
  fclose(fp); // close the file
}


// phase ii utility functions


void send_random_message_template(u8 M[HASH_INPUT_SIZE])
{ /* Send M immediately and clear the memory at the end  */

  MPI_Request requests[NSERVERS]; /* they will be used only here  */
  MPI_Status statuses[NSERVERS];
  // how about collective communications? I can't find a simple way using them.
  for (int i=0; i<NSERVERS; ++i) {
    MPI_Isend(M, HASH_INPUT_SIZE, MPI_UNSIGNED_CHAR, 0,
	      TAG_RANDOM_MESSAGE, MPI_COMM_WORLD, &requests[i]);
  }
  MPI_Waitall(NSERVERS, requests, statuses);
} /* clear stack variables */



void sender_process_task(u8 M[HASH_INPUT_SIZE],
			 WORD_TYPE Mstate[NWORDS_STATE])
{
  /* generate hashes and send them to a server as soon its buffer is complete */
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}
  // dgst := N bytes of the digest
  // ctr  := (this is a shortcut that allows us not to send the whole message)
  //---------------------------------------------------------------------------+

  int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                    + sizeof(CTR_TYPE); /* |dgst| + |ctr| */
  int server_number;


  // ---- Part1 : initializing sending buffers 
  MPI_Request request;
  /* int one_pair_size = sizeof(u8)*(N-DEFINED_BYTES) */
  /*                      + sizeof(CTR_TYPE); /\* |dgst| + |ctr| *\/ */

  u8* snd_buf = (u8*) malloc(one_pair_size
			     *PROCESS_QUOTA
			     *NSERVERS);

  /* Init attaching buffers for the buffered send */
  int buf_attached_size = NSERVERS*(MPI_BSEND_OVERHEAD + one_pair_size);
  u8* buf_attached = (u8*) malloc(buf_attached_size);
    

  /* Decide where to place the kth digest in server i buffer */
  size_t offset = 0; 
  const u64 ones = (1LL<<DIFFICULTY) - 1;

  /* pos i: How many messages we've generated to be sent to server i? */
  u32* servers_ctr = (u32*) malloc(sizeof(u32)*NSERVERS);
  /* set number of generated messages before sending to 0 */
  memset(servers_ctr, 0, sizeof(u32)*NSERVERS); 

  //-------------------------------------------------------------------------+
  //------- Part 2 : Generate hashes and send them  -------
  //-------------------------------------------------------------------------+
  MPI_Buffer_attach(buf_attached, buf_attached_size);
    
  while(1) { /* when do we break? never! */
    /* Find a message that produces distinguished point */
    find_hash_distinguished( M, Mstate, ones );

    //+ decide to which server to add to? 
    server_number = to_which_server((u8*) Mstate);



    /* 1st term: go to server booked memory, 2nd: location of 1st free place*/
    offset = server_number*one_pair_size + servers_ctr[server_number];
      
    /* save N-DEFINED_BYTES of MState in: snd_buf_dgst[offset] */      
    memcpy( ((u8*)Mstate) + DEFINED_BYTES, /* skip defined bytes */
	    (snd_buf+offset),
	    N-DEFINED_BYTES );

    /* record the counter, above we've recorded N-DEFINED_BYTES  */
    memcpy(M, 
	   ( (snd_buf+offset) + (N-DEFINED_BYTES) ),
	   sizeof(CTR_TYPE) );

      

    /* this server has one more digest */
    ++servers_ctr[server_number];

    if (servers_ctr[server_number] == PROCESS_QUOTA){
      /* we have enough messages to send to server (server_number) */
      MPI_Bsend_init(snd_buf,
		     PROCESS_QUOTA*one_pair_size,
		     MPI_UNSIGNED_CHAR,
		     server_number,
		     TAG_SND_DGST,
		     MPI_COMM_WORLD,
		     &request);
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
  // candidates. Save all candidates in a file called ?????                    |
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

  //---------------------------------------------------------------------------+
  // --- part 2: receive the initial inputs from all generating processors 
  //---------------------------------------------------------------------------+

  // process : 0 -> NSERVERS - 1 (receivers),
  //         : NSERVERS (archive)
  //         : NSERVERS + 1 -> nproc (senders)
  for (int i=0; i<nproc_snd; ++i) 
    MPI_Recv(initial_inputs+i*HASH_INPUT_SIZE,
	     HASH_INPUT_SIZE,
	     MPI_UNSIGNED_CHAR,
	     i + NSERVERS + 1,
	     TAG_RANDOM_MESSAGE,
	     MPI_COMM_WORLD,
	     statuses);


  //---------------------------------------------------------------------------+
  // --- Part 3: Receive digests and probe them
  //---------------------------------------------------------------------------+
  // enclose what's below within while loop to receive multiple times
  // Warning: not correct! nfound_collisions should be shared between all
  while (NNEEDED_CND_THIS_SERVER  > nfound_cnd) {
    //+ receive messages from different processors
  
    for (int i = 0; i<nproc_snd; ++i) {
      MPI_Irecv(rcv_buf+i*rcv_array_size,
		rcv_array_size,
		MPI_UNSIGNED_CHAR,
		i + NSERVERS + 1,
		TAG_SND_DGST,
		MPI_COMM_WORLD, requests);
    }

    //+ probe these messages 
    while (!flag) {/*receive from any server and immediately treat their dgst*/
      MPI_Waitsome(nproc_snd, requests, &outcount, indices, statuses);
      MPI_Testall(nproc_snd, requests, &flag, statuses);
      for (int j = 0; j<outcount; ++j){
	idx = indices[j]; /* number of sending process */
	buf_idx = idx*one_pair_size*PROCESS_QUOTA;
	msg_idx = idx*HASH_INPUT_SIZE;
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


void phase_ii(dict* d,
	      FILE* fp_possible_collisions,
	      size_t needed_collisions,
	      size_t difficulty_level)
{ 
  // some extra arguments to communicate with other
  
  // ==========================================================================+
  // Summary: Hash many different messages till we found enough collisions     |
  // --------------------------------------------------------------------------+
  // INPUTS:                                                                   |
  // `*d` : dictionary that has hashes from phase i                            |
  // `*fp_possible_collisions`: store all heashes h that dictionary tells they |
  //                            might exists. Also, store the corresponding    |
  //                            message.                                       |
  // NOTE: (64bytes message || 32bytes hash) (64bytes message || 32bytes hash) |
  // `needed_collisions` : How many messages we need that return positive hits |
  //                       when probing the dictionary with their digests.     |
  //---------------------------------------------------------------------------+
  // TODO:                                                                     |
  // - extra arguments to deal with other servers                              |
  // ==========================================================================+
  // process : 0 -> NSERVERS - 1 (receivers),
  //         : NSERVERS (archive)
  //         : NSERVERS + 1 -> nproc (senders)


  // --------------------- INIT MPI & Shared Variables ------------------------|
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
     
    /* init state */
    WORD_TYPE Mstate[NWORDS_STATE] = {HASH_INIT_STATE};

    /* set up a random message */
    /* M = 64bit ctr || 64bit nonce || random value */
    u8 M[HASH_INPUT_SIZE]; /* random word */
    /* Get a random message only once */
    CTR_TYPE* ctr_pt = (CTR_TYPE*) M; /* counter pointer  */
    getrandom(M, HASH_INPUT_SIZE, 1);
    ctr_pt[0] = 0; /* zeroing the first 64bits of M */

    
    /* Send the initial input to all receiving servers */
    send_random_message_template(M); 

    /* generate hashes and send them to a servers */
    sender_process_task(M, Mstate); /* never ends :) */
  }

  //+ todo receive processors
  //+ load dgsts from file to dictionary

  while (myrank < NSERVERS){ /* receiver, repeat infinitely  */ 
    //-------------------------------------------------------------------------+
    // I'm a receiving process: receive hashes, probe them, and send candidates 
    // Process Numbers: [0,  NSERVERS - 1]
    //-------------------------------------------------------------------------+


    
    /* Listen to senders till we accumulated enough candidates */
    receiver_process_task(d, myrank, nproc, nproc_snd);

    /* send the messages to the archive! */
    send_candidates_to_archive(myrank);
  }

  while (myrank == NSERVERS) /* I'm THE archive*/
    //+ todo truncate the messages file if it not multiple to HASH_INPUT_SIZE
    //+ this happen if the server was shut down during writing
    archive_receive();
      


  // The end that will never be reached in any case!
  MPI_Finalize();
  return;
}


// ------------------- Auxililary functions phase iii --------------------------
//+ todo complete these functions
size_t get_file_size(FILE *fp){
  /* return file size in bytes */
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  rewind(fp);
  return size;
} 


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

// phase iii
//+ master server has all the potential collisions messages and their digest
//+ a python script will combine them into a single file
//+ hash all message candidates
//+ master server sort the (hash, message) according to hash
//+ master server look at the long messag file, and search for
//+ hashes in the sorted list above





void phase_iii() /* single machine */
{
  // 1- load message candidates and hahs them
  // -- copy the hashes, and sort the copied hashes 
  // 2- load the middle states, hash the long message
  // -- after hashing a message block, query if the state is among the hashes
  //    from the candidates.
  
  
  // ----------------------------- PART 1 ------------------------------------
  // Pickle messages||digest in one array.
  /* load messages candidates, hash them, sort them */
  FILE* fp = fopen("data/receive/messages/archive", "r");
  size_t nmsgs = get_file_size(fp)/HASH_INPUT_SIZE;
  const WORD_TYPE state_init[NWORDS_STATE] = {HASH_INIT_STATE};
  
  // {(msg||dgst), ...} msg||dgst means concatenation
  u8* msgs  = (u8*) malloc( sizeof(u8)*nmsgs*HASH_INPUT_SIZE );
  u8* dgsts = (u8*) malloc(sizeof(u8)*nmsgs*(N) );
  u8* dgsts_orderd = (u8*) malloc(sizeof(u8)*nmsgs*(N) );
  
  fread(msgs, nmsgs, HASH_INPUT_SIZE, fp);
  fclose(fp);

  /* size_t step_size = HASH_INPUT_SIZE + N; */
  /* u8* msg_dgst = (u8*) malloc(nmsgs*step_size); */

  /* one thread is enough, we'll parellize it if it's a bottle neck */
  WORD_TYPE state[NWORDS_STATE];
  for (size_t i=0; i<nmsgs; ++i) {
    /* copy init state, then hash  */
    memcpy(state, state_init, NWORDS_STATE*WORD_SIZE);
    hash_single(state, &msgs[i*HASH_INPUT_SIZE]);

     /* get dgst in dgst */
    memcpy(&dgsts[i*N], state, N);
  }


  // ----------------------------- PART 3 ------------------------------------
  // sort msg_dgst according to the digests
  memcpy(dgsts_orderd, dgsts, nmsgs*N);
  qsort( dgsts_orderd, nmsgs, N, cmp_dgst);
  
  // ----------------------------- PART 2 ------------------------------------
  // hash the long message with each hashing probe 
  fp = fopen("long_message_middle_states", "r");
  size_t nmiddle_states  = get_file_size(fp)/(NWORDS_STATE*WORD_SIZE);
  /* How many state does a thread handle */
  size_t thread_load = nmiddle_states/omp_get_max_threads();
  

 #pragma omp parallel
  {
    void* ptr; // result of binary search
    int thread_num = omp_get_thread_num();
    size_t start = thread_load*thread_num;
    size_t end = thread_load*(1+thread_num);
    if (thread_num == omp_get_max_threads() - 1)
      end = nmiddle_states -1; /* last thread gets chunk + nmsg % nthreads */
    //+ todo how to incorporate counter?
    
    WORD_TYPE state_priv[NWORDS_STATE] = {HASH_INIT_STATE};
    u8 msg_priv[HASH_INPUT_SIZE] = {0};
    
    for (size_t i = start; i<end; ++i) {
      memcpy(state_priv, state_init, NWORDS_STATE*WORD_SIZE);
      // use counter to increment msg_priv
      hash_single(state_priv, msg_priv);
      ptr = bsearch(state_priv, dgsts_orderd, N, nmsgs, cmp_dgst);

      if (ptr) {
        /* remeber digest if the firs N bytes of state  */
	#pragma omp critical
	{
	
	size_t msg_idx = linear_search((u8*) state_priv, dgsts, nmsgs, N);
	u8* ptr_msg_collide = &msgs[msg_idx];
	

	
        printf("found a collision at %lu\n", i);
	printf("hash long message at %lu:\n", i);
	print_byte_array((u8*) state_priv, N);
	
	printf("while the the following message:\n");
	print_byte_array(ptr_msg_collide, HASH_INPUT_SIZE);
	puts("produce the following hash:");
	print_byte_array(&dgsts[msg_idx], N);
	
	}// end critical region
      } // end if condition
      
    
    } // end for loop, thread's main work

  } // end parallel region
} // quit the function
