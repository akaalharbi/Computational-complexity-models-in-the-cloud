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


void send_random_message_template(u8 M[HASH_INPUT_SIZE])
{ /* Send M immediately and clear the memory at the end  */
  
  MPI_Request requests[NSERVERS]; /* they will be used only here  */
  MPI_Status statuses[NSERVERS];


  // how about collective communications?
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
  // I am a sending processor. I only generate hashes and send them.         |
  // Process Numbers: [NSERVERS + 1,  nproc]                                 |
  //-------------------------------------------------------------------------+

  // ----------------------------- PART 0 --------------------------------- //
  // 1- set up initial random message.
  // 2- set up counting pointer.


  /* M = 64bit ctr || 64bit nonce || random value */
  u8 Ms[16][HASH_INPUT_SIZE]; /* random word */
  u8 M[HASH_INPUT_SIZE]; /* random word */
  
  /* Get a random message only once */
  CTR_TYPE* msg_ctr_pt = (CTR_TYPE*) M; /* counter pointer  */
  getrandom(M, HASH_INPUT_SIZE, 1);
  msg_ctr_pt[0] = 0; /* zeroing the first 64bits of M */

  for (int i = 0; i<16; ++i) {
    memcpy(Ms[i], M, HASH_INPUT_SIZE);
  }
  
  char txt[50];
  snprintf(txt, sizeof(txt), "sender #%d template=", myrank);
  print_byte_txt(txt, M,HASH_INPUT_SIZE);

  // ----------------------------- PART 1 --------------------------------- //
  // 1-  Sen the initial input to all receiving servers 
  send_random_message_template(M); /* send template to all servers */



  // ------------------------------ PART 2 ----------------------------------- +
  // 1- generate disitingusihed hashes                                         |
  // 2- decide which server should probe the disitingusihed hash               |
  // 3- if server i buffer has enough hashes, send them immediately.           |
  //---------------------------------------------------------------------------+
  // snd_buf ={ dgst1||ctr1, ..., dgst_k||ctr_k}                               |
  // dgst := N bytes of the digest                                             |
  // ctr  := (this is a shortcut that allows us not to send the whole message) |
  //---------------------------------------------------------------------------+
  WORD_TYPE Mstate[NWORDS_STATE] = {HASH_INIT_STATE};
  size_t one_pair_size = sizeof(u8)*(N-DEFINED_BYTES)
                       + sizeof(CTR_TYPE); /* |dgst| + |ctr| - |known bits|*/

  int server_number = -1;
  // { (server0 paris) | (server1 pairs) | ... | (serverK pairs) }
  u8* snd_buf = (u8*) malloc(one_pair_size
			     *PROCESS_QUOTA
			     *NSERVERS);

  /* Decide where to place the kth digest in server i buffer */
  // How many bytes are reserverd for each server in snd_buf
  size_t nbytes_per_server = one_pair_size * PROCESS_QUOTA;
  size_t offset = 0; /* which index within a  server buffer should we pick */

  // use the mask to decide 
  const u64 mask_test = (1LL<<DIFFICULTY) - 1; /* mask_test & diges =?= 0 */

  /* pos i: How many messages we've generated to be sent to server i? */
  u64 servers_ctr[NSERVERS] = {0};
  int snd_ctr = 0; /* how many messages have been sent */


  
  // ------------------------------ PART 2 ----------------------------------- +
  // Generate hashes and send them 
  double time_start = wtime(); 
  printf("sender #%d: done init mpi, and sharing the its template."
	 "going to generate hashes \n", myrank);
  
  // find_hash_distinguished_init(); 



  while(1) { /* when do we break? never! */
  /* while(i<1) { /\* when do we break? never! *\/ */

    /* Find a message that produces distinguished point */
    find_hash_distinguished(Ms, M, Mstate, msg_ctr_pt, mask_test);

    //+ decide to which server to add to? 
    server_number = to_which_server((u8*) Mstate);


    /* 1st term: go to server booked memory, 2nd: location of 1st free place*/
    offset = server_number * nbytes_per_server
           + servers_ctr[server_number] * one_pair_size;
    // recall que one_pair_size =  |dgst| + |ctr| - |known bits|




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
      

    
    if (servers_ctr[server_number] >= PROCESS_QUOTA){

      printf("===============================================\n"
             "sender #%d -> recv #%d before sending %0.4fsec\n"
	     "===============================================\n\n",
	     myrank, server_number, wtime() -  time_start );

      time_start = wtime();
      MPI_Send(&snd_buf[server_number*nbytes_per_server],
		PROCESS_QUOTA*one_pair_size,
		MPI_UNSIGNED_CHAR,
		server_number,
		TAG_SND_DGST,
		MPI_COMM_WORLD);

      
      ++snd_ctr;

      printf("-----------------------------------------------\n"
             "sender #%d -> recv #%d sending done %0.4fsec\n"
	     "-----------------------------------------------\n\n",
	     myrank, server_number, wtime() -  time_start );


      // debugging: check if sender and receiver have the same digest
      /* if (server_number == 7 && myrank == 9) { */
      /* 	snprintf(txt, sizeof(txt), "sender#%d,  server=%d, snd_buf=", */
      /* 		 myrank, server_number ); */
      /* 	print_byte_txt(txt, */
      /* 		       &snd_buf[server_number*nbytes_per_server], */
      /* 		       one_pair_size*PROCESS_QUOTA); */
      /* } */


      /* time_start = wtime(); */

      /* It is enough to reset the counter. The memroy will be rewritten */  
      servers_ctr[server_number] = 0;
      time_start = wtime();
    }
  } 


  free(snd_buf);
  return; // au revoir.

} // MPI_Finalize
