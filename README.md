# Long message second pre-image attack on SHA2-256 truncated to at most 96-bits


## Overview

An implementation of long message attack on a truncated SHA2-256 in parallel using message passing interface (MPI) and AVX512, the user can choose the truncation amount.


The goal of this project is to estimate the cost of the attack according to various benchmarks( cpu hours, energy, ...) and comparing them with the theoritical cost model. Based in our experience with [Grid5000](https://grid5000.fr), the booked RAM expands with the number of requested node however this doesn't apply to hard-disk

### Attack description 

 
This is a generic attack on any hash function based on Merkle-Damgard construction without strengthening with a compression function $f:\{0, 1\}^{m+k} -> \{0, 1\}^n$

If we have `M = M0 | M1 | \dots | M_L ` where `|M_i| = m `, i.e. the message M is broken into blocks of the same length, $m$ bits. Then, $h(M) = S_l$ where $S_i = f(M_i, S_{i-1})$ where $S_{-1} = IV$  for a fixed $IV$.


**Warning** Since it's a proof of concept, we got lazy and trucated only the last state. If we truncate all the middle states while hashing the long message, i.e. edit `hash_single` and wrap the `sha256_multiple_x16_tr` to truncate their output states then the program would work as the same.


We hash a long enough message $M$ and store all the middle states $S_i$, denote the number of blocks in $M$ as $L = 2^{l}$. Then, hash a random message $h(R) = sha2-256\left(R, IV\right) $. If $h(R) = S_i$ then `R| M_i | \dots | M_l` is a second preimage of $M$.

In our case the compression function: $sha2-256: \{0, 1\}^{512 + 256} -> \{0, 1\}^{256}$.

### Best Theoritical Attack Parameters
If the digest lenght is $n$-bits then choosing $l := 2^{\frac{n}{2}}$ minimizes the required number of expected random hashing to $2^{\frac{n}{2}}$

### Best Attack Parameters in the Real World
todo!


## How does the program work?


### phase i: hash a long message
The long message `M := M0 | ... | M_i | \dots ` is defined as `M_i := i | 0 \dots 0`, this avoids being trapped into repeated cycle if we always hash the same message although it's unlikely event.

Storing all middle states require a large hard-disk, more than a petabyte, thus we store the middle state after $2^{30}$ hashes in the file `data/states`. We can recover all the middle states in parallel.

This phase requires a cpu that support Intel SHA extensions.


### phase ii (under construction): 
This is the main working horse. It will reconstruct the long message in parallel, adding the middle states to a dictionary, then generates many random message and collects those that match partially an entry in the dictionary.

todo: what is stored as index and as a value in the dictionry, dictionary probing method, probabilities of false negative, false positive, cache misses, huge pages and tlb misses

todo: breaking processes into senders and receivers, sha2-256-16way, MPI message structure, overview of the main steps in phase ii, adjusting speed of senders, 

todo: estimating the needed number of positive probes, explain why an increase number of needed positive probes doesn't mean necessairly higher number of hashing, squeezing all memory bits


### phase iii (under construction)


## Getting started

### perquisites

- nasm assembler
- GCC (probably >= 7.3.1, tested on 10.2.1-6)
- CPU supports AVX-512 needed by `phase_ii`
- CPU has Intel SHA extensions for `phase_i` & `phase_iii` (the latter might require `avx512`

### Run
We assume all nodes have same amount of RAM. 

A typical run example:

```
python run.py  --nservers 8 --receivers 8 -N 12 --ram 64000000000 --interval 30 --difficulty 2
```

> todo: remove the lines (they are specific to [Grid5000.fr](https://www.grid5000.fr/w/Grid5000:Home)
```
os.system("echo 'core' | sudo-g5k tee /proc/sys/kernel/core_pattern")
os.system("sudo-g5k apt install nasm")
```
> todo: more instruction on how to define the mpi's machinefile (currently we're using $OAR_NODEFILE as machinefile)

- `nservers`: how many nodes will be used?
- `receivers`: The total number of receivers, there should be at least one receiver per node. 
- `N`: How many bytes of digest we are going to consider.
- `ram`: How much RAM in one node should be allocated for receivers in that node.
- `interval`: laissez tomber
- `difficulty`: How many bits in the digest should be zero. This mean the sender will only send hashes that have exactly `difficulty` zeros at the end.



todo: add benchmarks section
