# Long message attack
- Phase I  : Construct a long message
- Phase II : Find sufficient number of potentiall collisions. 
- Phase III: Filter all false positives, record messages the produces collisions.

# Known Bugs:
- Phase II: doesn't work

## todo (Techincally not bugs):
- Rewrite the format of the files
- restructure files in the folder
- add time stamp for the files
- compute how many potentiall collisions should be computed in phase ii
- write phase iii
- edit benchmark
- edit verifyo
- when writing the long message, we might record several internal states of it wiht index in a file. This allows parallel execution in phase iii
- `log/` folder is for statistics and testing output
- `messages/` there should be a distinction between potential messages and messages the leads to collisions

# Performance Evolution
In these measurements, we used a predictable prng in order to test the same numbers on all variants. Also, all parameters below have no cycle in Phase I. We believe this is realistic model.

## n=52, l=25
- 13 oct 2022:  13.16 sec (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
- 27 sep 2022: 976.07 sec (ɑ=0.50, sha256, openmp)

## n=51, l=25
- 13 oct 2022:   7.41 sec (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
- 27 sep 2022: 322.98 sec (ɑ=0.50, sha256, openmp)

## n=50, l=25
- 13 oct 2022:  4.22 sec (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
- 27 sep 2022: 51.15 sec (ɑ=0.50, sha256, openmp)


## n=44, l=18
- 13 oct 2022:   3.20 sec (ɑ=0.90, sha256-x86, simd, openmp, 64bit key store max)
- 27 sep 2022: 103.01 sec (ɑ=0.50, sha256, openmp)


# Summary:
This folder has three main executables `long_message_attack`, `verify_hash', `bench_mark` which generate two colliding messages, verifies the results, and benchmark the dictionary and sha256 respectively. 

# definitions
long message: a large message of zeros.

## Usage:
Basic usage:
make clean && make

./long_message_attack n l

This will mount an attack of truncated sha256 to n-bits using a long message of length 2^l. The random message that collides with the long message can be found in `messages/n_l` in binary format (read it as unsigned char[64]). Also, you can find benchmark and the index where the collision occurs at `statistics_parallel/n_l.txt`.  Substitute n and l, if l>= n/2 then there is a high chance of having a cycle. In this case, you can only find a file in `statistics_parallel/n_l.txt`


### Verifies the results
./verify_hash "messages/n_l" n
This will read the file n_l, truncate the output to n bits, and search for an index where the message found in `messages/n_l` collides with the long message. If the index mentioned in "messages/n_l" doesn't produce a collision it will indicate that they are not colliding.




# compile:
run `make` and the executable `long_message_attack` is the a parallel version of the attack, apologies for the confusion with what written below.

## Flags
`-fopenmp` for parallel executable (only Phase II can be parallel)
`-DVERBOSE_LEVEL` print intermediate information (it will flood the program use with care)


- parallel
gcc -fopenmp -O3 long_message_attack.c sha256.c dict.c util/util_char_arrays.c util/memory.c -lm -o long_message_attack_parallel

// DEBUG
gcc -DVERBOSE_LEVEL -fopenmp -O3  -g long_message_attack.c sha256.c dict.c util/util_char_arrays.c util/memory.c -lm -o long_message_attack_parallel




- serial
gcc -O3  -g long_message_attack.c sha256.c dict.c util/util_char_arrays.c util/memory.c -lm -o long_message_attack

// with debug using
gcc -DVERBOSE_LEVEL -O3  -g long_message_attack.c sha256.c dict.c util/util_char_arrays.c util/memory.c -lm -o long_message_attack



- verify hash
gcc -g verify_hash.c sha256.c -lm -o verify_hash

cost statistics:
sudo perf stat -a  -e power/energy-pkg/ ./long_message_attack_parallel 63 29

sudo perf stat -o stats.txt -a -r 2 -e power/energy-pkg/ ./long_message_attack 30 15


sudo perf stat -a -r 2 -e power/energy-pkg/ ./long_message_attack 30 15

usage:
./longa_message_attack n l
n: 0<n<257 which is the number of bits in the output of the compression function
l:float 2^l is the number of blocks

- the two messages, that collide, will be saved in the files `message1` `message2`


files:
- 
- statistics*/n_l_bits.txt files contain the system power usage while the program was running.
- long_message_attack.c: main file that implements the attack
- verify_hash.c: reads 64 bytes from a binary file then it hashes a long message till a collision is found.
                 returns the index idx where the collision occurs.
- dict.*: chained dictionary implementation using linked list
- util_char_arrays.* : some useful functionalities
- ignore the tmp file



taken from the internet: 
- sha256.* (modified)
- crypto-algorithms-master.zip
- memory.c

questions:


issues:


future improvements:


possible optimizations:
- profile the code


Benchmarks:
- Number of sha256 transforms ~2^21 calls / sec
- Memory writes ~100MB/sec
