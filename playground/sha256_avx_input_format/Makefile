lib:
	cd sha256_intel_avx && make clean && make all
	


main: main.c
	gcc  main.c -L./sha256_intel_avx -L./sha256_ni \
	-lsha256_avx -lsha256_ni -o main


all: main lib

.PHONY: clean lib

clean:
	rm -f *.o *.a main


