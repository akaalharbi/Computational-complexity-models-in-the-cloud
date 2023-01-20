NASM_FLAGS := -Werror -felf64 -Xgnu -gdwarf -DLINUX -D__linux__ $(NASM_INCLUDES)
CFLAGS = -Wall
LDLIBS = -lm


timing.o: timing.c

vsha256.o: vsha256.c
	gcc -c -g -march=native vsha256.c -o vsha256.o	

driver.o: arch_avx2_type1.h driver.c

sha256_oct_avx2.o: sha256_oct_avx2.asm
	nasm -o $@  $(NASM_FLAGS) $<

sha256_ni.o: sha256_ni.c
	gcc -c -msse4.1 -msha sha256_ni.c -o sha256_ni.o

main: driver.o sha256_oct_avx2.o sha256_ni.o timing.o vsha256.o

all: main

.PHONY: clean

clean:
	rm -f *.o main
