NASM_FLAGS := -Werror -felf64 -Xgnu -gdwarf -DLINUX -D__linux__ $(NASM_INCLUDES)
CFLAGS = -Wall


# driver: driver.o sha256_oct_avx2.o
# avx2
c_sha256_oct_avx2.o: arch_avx2_type1.h c_sha256_oct_avx2.c

sha256_oct_avx2.o: sha256_oct_avx2.asm
	nasm -o $@ $(NASM_FLAGS) $<




# avx512
c_sha256_x16_avx512.o: arch_avx512_type1.h c_sha256_x16_avx512.c

sha256_x16_avx512.o: sha256_x16_avx512.asm
	nasm -o $@ $(NASM_FLAGS) $<




# final library
libsha256_avx.a: c_sha256_x16_avx512.o sha256_x16_avx512.o c_sha256_oct_avx2.o sha256_oct_avx2.o
	ar -rcv $@ $^






all: libsha256_avx.a

.PHONY: clean

clean:
	rm -f *.o *.a
