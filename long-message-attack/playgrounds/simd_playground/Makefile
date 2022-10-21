LDLIBS  = -lm
CFLAGS  = -g -march=native
objects = simd_play.o



all: simd_play

simd_play: $(objects)

# objects and their dependencies


.PHONY: clean
clean:
	rm -f $(objects)
	rm -f simd_play

