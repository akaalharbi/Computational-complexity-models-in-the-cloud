CFLAGS = -Wall
LDLIBS = -lm


timing.o: timing.c

vsha256.o: vsha256.c
	gcc -c -g -march=native vsha256.c -o vsha256.o	


lib:
	cd sha256_intel_avx && make clean && make all

main: timing.o vsha256.o
	gcc  main.c -L./sha256_intel_avx -L./sha256_ni \
	-lsha256_avx -lsha256_ni -lm  timing.o vsha256.o -o main


all: main lib

.PHONY: clean

clean:
	rm -f *.o main
