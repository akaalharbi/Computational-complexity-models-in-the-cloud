#CC=clang
#CC=gcc
LDLIBS  = -lm
INCLUDE = include
INC = -I$(INCLUDE)

ifeq ($(CC), clang)
	CFLAGS = -g -O3 -fopenmp -Wall -march=native -msha  #-Rpass-analysis=loop-vectorize -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Xanalyzer -analyzer-constraints=z3
else
	CFLAGS  = -g -O3  -march=native -std=c11  #-DVERBOSE_LEVEL=2 
endif



# store *.o files in obj/
OBJDIR = obj

# all source files can be found in these folders
SRC = .
SRC_UTIL = $(SRC)/util
#SRC_SHA256 = $(SRC)/sha256

## extract all *.c filenames from the directories
FILENAMES  =  $(wildcard $(SRC)/*.c)
FILENAMES +=  $(wildcard $(SRC_UTIL)/*.c)

# Substitution References: replaces all *.c with *.o
# note that it will keep the directory structure
OBJECTS := $(FILENAMES:$(SRC)/%.c=$(OBJDIR)/%.o)

# save objects in 
#OBJECTS := $(addprefix $(OBJDIR)/, $(OBJECTS))
# $(info FILENAMES = $(FILENAMES))
# $(info OBJECTS = $(OBJECTS))
# $(info CC=$(CC))
# # list of all *.o


# BUILD OBJECT FILES IN OBJECTDIR directory
$(OBJDIR)/%.o: $(SRC)/%.c
	mkdir -p '$(@D)'
	$(CC) -c $< $(INC) $(CFLAGS)  -o $@






# todo make complete creating two targets


all: main
	mkdir -p obj



# remove all $(TARGETS) members from dependencies
# and add long_message_attack.o as a dependency
main: $(OBJECTS)
	$(CC)  $? $(LDFLAGS) $(LDLIBS) -o $@ 


.PHONY: clean
clean:
	rm -f $(OBJECTS)

