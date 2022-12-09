// dummy file that will be removed later
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
