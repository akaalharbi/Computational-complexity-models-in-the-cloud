NASM_FLAGS := -Werror -felf64 -Xgnu -gdwarf -DLINUX -D__linux__ $(NASM_INCLUDES)
CFLAGS = -Wall
LDLIBS = -lm





sha256_ni.o: sha256_ni.c
	gcc -c -msse4.1 -msha sha256_ni.c -o sha256_ni.o



libsha256_ni.a: sha256_ni.o
	ar -rcv $@ $^

all: libsha256_ni.a

.PHONY: clean

clean:
	rm -f *.o *.a
