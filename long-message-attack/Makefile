#CC=g++
LDLIBS  = -lm
CFLAGS  = -g -fopenmp -Wall -march=native -msse4 -msha -std=c11 -fopt-info-vec -fopt-info-omp-vec-optimized-missed  #-DVERBOSE_LEVEL=2 
LDFLAGS = -fopenmp -march=native -msse4 -msha
objects = long_message_attack.o sha256-x86.o dict.o util/memory.o util/util_char_arrays.o



all: long_message_attack

long_message_attack: $(objects)

# objects and their dependencies
long_message_attack.o: sha256-x86.h types.h util/util_char_arrays.h shared.h
sha256-x86.o: sha256-x86.h
dict.o: dict.h types.h shared.h
util/memory.o: util/memory.h
util_char_arrays.o: util/util_char_arrays.h


.PHONY: clean
clean:
	rm -f $(objects)
	rm -f long_message_attack
