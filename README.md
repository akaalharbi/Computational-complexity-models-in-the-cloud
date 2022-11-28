# Long message attack
This repository is still a work in progress that should implement long message attack on parallel using multiple computers. The goal is to attack 96-bits (or higher!) of `sha-256` output. 

## Attack Overview
For a better description see (find somone explains lma)[]. The attack can be parametrized by $n, l$ where denote $n=nbits$ to be attacked and $2^l$ is the number of elements stored in `Phase I` (number of blocks in the long message)


- Phase I  : Constructs a long message
- Phase II : Finds sufficient number of potentiall collisions. 
- Phase III: Filter all false positives, record messages the produces collisions.

## Project milestones
1. Mounting the attack on a single core on local machine, find highest $n$ and $l$
2. Mounting the attack on a single using all cores, find highest $n$ and $l$
  - Estimate the resources needed to attack $n:=96$, in our case it seemed to be 1024 servers!
3. Optimize the attack on single machine. Use the largest possible memory, and the fastest sha256 evaluation.
   - Estimate again the needed number of servers (assuming that all servers are like the tall maching cluster.lip6.fr)
4. Mounting the attack on multiple devices 


We are in ** step 3.**


## Message structure
sha256 has a 
uint32_t state[8] =  {A, B, C, D, E, F, G, H};


Casting state to uint64_t:

(uint64_t*) state; //= { (B<<32|A), (D<<32|C), (F<<32|E), (H<<32)|G } 


uint32_t digest[3] = {A, B, C}; // 96 bits

### Distinguished point and server number
A = rest||server_number||bits for distinguished point

### What to store in the dictionary

the index is 
Let m:= nbits dedicated for distinguished points + nbits dedicated for server_number

idx= (B<<32|A) >> m // remove truncate distinguished points and server number
idx = idx mod (nbuckets)
Store whole word C (or state[2]) as a value in the dictionary 





# Usage


`make clean && make`

Then run it
`./long_message_attack n l` where $n$ and $l$ are the attack parameters.



# Work to be done:

## bugs: 
- 

## todo (Techincally not bugs):
- change the counter in `void find_hash_distinguished` to 128bits
  i.e. In phase II, start with random message then increment 128 bit by one each time
- move hard coded choices to config.h, e.g. nbits to be stored as an index
- Due to optimizations introduced, the current code is not flexible with the choice of n!
-- Restore the flexibility, so that we can test it on small n.
- Stick to 80 characters per line. 
- Draw communication model




# Performance Evolution
In these measurements, we used a predictable prng in order to test the same numbers on all variants. Also, all parameters below have no cycle in Phase I. We believe this is realistic model. These tests have been conducted on the same machine with the same OS.

## sha256 performance:
- sha256-x86 elapsed=1.046995sec, 32048322.000000 hash/sec≈2^24.933745 
- parallel sha256-x86 elapsed=1.737093sec, 386328576.000000 hash/sec≈2^28.525253 
- vsha256 elapsed=1.057237sec, 31737852.000000 hash/sec≈2^24.919701 
- vsha256 parallel eapsed=3.001098sec, 223614368.000000 hash/sec≈2^27.736438 


## to gather or not to gather:
Benchmarking dictionary lookup where 4-elements simultaneously using gather instruction, the look up are slower than linear lookup of a single element. On the other hand, when benchmarking `long\_message\_attack 52 25` the gather version seems to be 20% faster. 

These conflicting measures can be resolved by observing the two facts:
1. The old code has a dictionary size almost 2x larger than dictionary used in gather version. This is because, we decided to remove values array. 
2. Counting the number of trials, due to to the prng, gather version was lucky when it found the collision after ~30M trials. The linear version had to go through ~70M trials and was around ~20% slower!

## Verdict:
Gather instruction is not worth using in our case. 

### Addendum:

I found claims that:
	- Gather is 8x faster on new amd threadripper (see https://twitter.com/ChrisGr93091552/status/1562583460992397317 )
	- Gather can be as fast as linear scan (see to use or not to use the simd gather https://dl.acm.org/doi/abs/10.1145/3533737.3535089 )


## Long message attack on a single machine:
### n=52, l=25
- 13 oct 2022(updated):  5,9917 +- 0,0192 seconds time elapsed  ( +-  0,32%) (ɑ=0.90, sha256-x86, simd(single element search), openmp, 64bit key store)
- (not-accurate) 27 sep 2022: 976.07 sec (ɑ=0.50, sha256, openmp)
- (Warning: This result is faster due to incorrect behavior! when 3/4 decides to end the probe, ) 03 nov 2022: 3,9475 +- 0,0131 sec (ɑ=0.9, sha256-x86, simd(multiple search, stop if 3/4 found empty slot), openmp, 64bit key store) 




### n=51, l=25

- 02 nov 2022:  5,0265 +- 0,0145 sec (ɑ=0.9, sha256-x86, simd(multiple elements search), openmp, 64bit key store)
- 13 oct 2022(updated):  5,611 +- 0,0292 sec  (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
- 27 sep 2022(not-accurate): 322.98 sec (ɑ=0.50, sha256, openmp)
- (Warning: This result is faster due to incorrect behavior!) 03 nov 2022: 3,9668 +- 0,0143 sec (ɑ=0.9, sha256-x86, simd(multiple search, stop if 3/4 found empty slot), openmp, 64bit key store) 
### n=50, l=25
(c rand_r, seed is the thread number)
- 02 nov 2022: 4,0711 +- 0,0105 secc (ɑ=0.9, sha256-x86, simd(multiple elements search), openmp, 64bit key store)
- 13 oct 2022:   5,6730 +- 0,0871 sec (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
(random numbers)
- (Warning: This result is faster due to incorrect behavior!) 03 nov 2022: 3,95056 +- 0,00790 sec (ɑ=0.9, sha256-x86, simd(multiple search, stop if 3/4 found empty slot), openmp, 64bit key store) 
- 27 sep 2022(not-accurate): 51.15 sec (ɑ=0.50, sha256, openmp)


### n=44, l=18
(todo restore old data)
- 03 nov 2022: 0,01836 +- 0,00231 sec (ɑ=0.9, sha256-x86, simd(multiple search, stop if 3/4 found empty slot), openmp, 64bit key store) 



