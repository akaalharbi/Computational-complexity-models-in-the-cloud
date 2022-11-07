LDLIBS  = -lm
CFLAGS  = -g -O3 -Wall -march=native -msse4 -msha -fopenmp
LDFLAGS = -msse4 -msha -fopenmp
objects = sha256wrapper.o sha256_one_block_avx.o sha256-x86.o #sha256_mult_avx.o



all: sha256wrapper
sha256wrapper: $(objects)
sha256-x86.o: sha256-x86.h

.PHONY: clean
clean:
	rm -f $(objects)
