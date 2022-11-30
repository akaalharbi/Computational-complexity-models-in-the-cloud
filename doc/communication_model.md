Warning: all variables in this document are local to this document. e.g.
we used m as the number of servers while in `README.md` m denoted sth else.

Let's say we have m computers/servers that communicate with each other directly.
1. Assumption: each thread within a server contributes equally to the generated hashes.
2. Assumption: each server sends a fixed number of messages each time.

The first assumption enables each thread to insert elements in send buffer
without using locks. The second assumption facilitates easier sectioning of
receive without explicit synchronizations between servers.






# High level model:


MPI_Ireceive <--
Generate enough messages
MPI_Wait() // for receiving messages 


MPI_Isend -->
probes the local dictionary with received messages
MPI_Wait() // for sending messages


Also, use `MPI_All_reduce` to update number of founded potential collisions.


# Messages Format


## INIT
each server sends its initial random message:
if we have m servers, then each server will have a copy of
- {rd_1, rd_2, rd_3, ..., rd_m}




## Receiver point of view:

server i will receive sth like: 
- {  hash_1,   hash_2, ...,     hash_n}
  {offset_1, offset_2, ....,  offset_n}

offset_i : 128bit
Receiving server can construct message related to hashi by setting the last
4 words of mj to offset where mj is the initial random message of server j.
j depends on i.

### Determining j, or what server has generated the message:
Say we are in server t, when it probes a hash_i and the dictionary tells this
digest exists we wish to record the message the generated digest h.

Recovering the message requires computing (set last 4 words of m\_j to offset\_j)
for unkown j yet.  Since each server writes on specified number of locations, 


e.g. if each server contributes equally we pick j s.t. i \in [j*n/m, (j+1)n/m )
or simply j = floor( i / (n/m) ). If servers don't contribute equally then

### Servers don't contribute equally
If there are m servers, server k sends r_k messages.

Then construct the sorted list:
R = {r_1, r_1+r_2, ..., r_1 + r_2 + ... + r_m }

to find j:
find the least index j where R[j] > i.




## Sender point of view:

The sender k will generate simultaneously r_k messages for each sever. 
The message generations will be done in parallel using server k parallel
available threads. For now, we will reserver the master thread to communicate
with other servers.





### Messages generations
The server k has to send r_k messages to each of m servers (including self). 
It has 
snd_buf_1[r_k]

snd_buf_2[r_k]

...

snd_buf_m[r_k]



snd_buf_i[r_k] is divided into equal number of chunks as :

snd_buf_i[r_k/nthreads][nthreads % r_k]

e.g. thread 0 will only write entries on snd_buf_i[0][_]

The two brackets notation is converted into:


snd_buffer_i[ t*(r_k/nthreads) + _ ] where t is the thread number.


i is chosen based of computed digest, see Message Structure in `README.md`



# Techincal details:

## size of send and receive buffers
This is still undetermined.




---


