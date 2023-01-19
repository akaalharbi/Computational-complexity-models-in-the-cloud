NASM_FLAGS := -Werror -felf64 -Xgnu -gdwarf -DLINUX -D__linux__ $(NASM_INCLUDES)
CFLAGS = -Wall
all: driver

driver: driver.o sha256_oct_avx2.o

driver.o: arch_avx2_type1.h driver.c

sha256_oct_avx2.o: sha256_oct_avx2.asm
	nasm -o $@ $(NASM_FLAGS) $<

.PHONY: clean

clean:
	rm -f *.o driver