NASM_FLAGS := -Werror -felf64 -Xgnu -gdwarf -DLINUX -D__linux__ $(NASM_INCLUDES)
CFLAGS = -Wall


# driver: driver.o sha256_oct_avx2.o

c_sha256_oct_avx2.o: arch_avx2_type1.h c_sha256_oct_avx2.c

sha256_oct_avx2.o: sha256_oct_avx2.asm
	nasm -o $@ $(NASM_FLAGS) $<


libsha256_oct_avx2.a: c_sha256_oct_avx2.o sha256_oct_avx2.o
	ar -rcv $@ $^


all: libsha256_oct_avx2.a

.PHONY: clean

clean:
	rm -f *.o *.a
