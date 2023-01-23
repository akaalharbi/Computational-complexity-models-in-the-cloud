


main: main.c
	gcc  main.c -L./sha256_intel_avx -L./sha256_ni -lsha256_oct_avx2 -lsha256_ni -o main


all: main

.PHONY: clean

clean:
	rm -f *.o *.a main
