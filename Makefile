##############################################
#             Table of Contents              #
# 1: Basic Flags in the compile and linking  #
# 2: where to put and to find files?         #
# 3: Compiling and Linking                   #
##############################################


# Part 1: Basic Flags in the compile and linking
CC=gcc

# Doesn't work
# MPI_INCLUDE := $(shell mpicc -showme:compile | grep -e . )
# MPI_LINK := $(shell mpicc -showme:link)


MPI_INCLUDE = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
MPI_LINK = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi



INCLUDE = -Iinclude $(MPI_INCLUDE)
LDLIBS  =  $(MPI_LINK) -lm -L./lib/sha256_intel_avx/ -lsha256_avx

LDFLAGS = -fopenmp -pthread -O3 -flto -g  # -fsanitize=address
# note: -fanalyzer doesn't report when -flto is enabled
CFLAGS =  -O3 -flto -g  -fopenmp -Wall -march=native -msha  -std=c11 #-fopt-info-all -g -fanalyzer -fsanitize=address 
#CFLAGS += -DVERBOSE_LEVEL=2 



# Part 2: where to put and to find files?
# store *.o files in obj/
OBJDIR = obj

# all source files can be found in these folders
SRC = src
SRC_UTIL = $(SRC)/util
#SRC_SHA256 = $(SRC)/sha256

## extract all *.c filenames from the directories
FILENAMES  =  $(wildcard $(SRC)/*.c)
FILENAMES +=  $(wildcard $(SRC_UTIL)/*.c)

# Substitution References: replaces all *.c with *.o
# note that it will keep the directory structure
OBJECTS := $(FILENAMES:$(SRC)/%.c=$(OBJDIR)/%.o)





# Part 3: Compiling and Linking


lib:
	cd lib/sha256_intel_avx && make clean && make all

# BUILD OBJECT FILES IN OBJECTDIR directory
$(OBJDIR)/%.o: $(SRC)/%.c
	mkdir -p '$(@D)'
	$(CC) -c $< $(INCLUDE) $(CFLAGS)  -o $@




TARGETS = test phase_i phase_ii phase_iii

# REMOVE TARGETS FROM $(OBJECTS)
TARGET_OBJECTS = $(addprefix $(OBJDIR)/,  $(addsuffix .o, $(TARGETS)) )
COMMON_OBJECTS = $(filter-out $(TARGET_OBJECTS), $(OBJECTS) )

$(info common objects = $(COMMON_OBJECTS))

all: $(TARGETS) lib
	mkdir -p obj
	mkdir -p data
	mkdir -p data/digests
	mkdir -p data/messages
	mkdir -p data/digests
	mkdir -p data/messages
	mkdir -p data/stats




# we wish to build X:
# 1- remove all $(TARGETS) members from dependencies
# 2- add X.o as a dependency

test: $(OBJDIR)/test.o $(COMMON_OBJECTS)
	$(CC)  $^ $(LDFLAGS) $(LDLIBS) -o $@ 


phase_i: $(OBJDIR)/phase_i.o $(COMMON_OBJECTS)
	$(CC)  $^ $(LDFLAGS) $(LDLIBS) -o $@ 

phase_ii: $(OBJDIR)/phase_ii.o $(COMMON_OBJECTS)
	$(CC)  $^ $(LDFLAGS) $(LDLIBS) -o $@ 

phase_iii: $(OBJDIR)/phase_iii.o $(COMMON_OBJECTS)
	$(CC)  $^ $(LDFLAGS) $(LDLIBS) -o $@ 


.PHONY: clean
clean:
	rm -f $(OBJECTS)
	rm -f $(TARGETS)
purge:
	rm -f $(OBJECTS)
	rm -f $(TARGETS)
	rm -rf data
