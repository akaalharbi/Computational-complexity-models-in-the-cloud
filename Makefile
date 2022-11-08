#CC=g++
LDLIBS  = -lm
CFLAGS  = -g -O3 -fopenmp -Wall -march=native -msse4 -msha -std=c11 -fanalyzer -fopt-info-vec -fopt-info-omp-vec-optimized-missed #-DVERBOSE_LEVEL=2 
LDFLAGS = -fopenmp
INCLUDE = include
INC = -I$(INCLUDE)

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

# save objects in 
#OBJECTS := $(addprefix $(OBJDIR)/, $(OBJECTS))
# $(info FILENAMES = $(FILENAMES))
# $(info OBJECTS = $(OBJECTS))

# # list of all *.o


# BUILD OBJECT FILES IN OBJECTDIR directory
$(OBJDIR)/%.o: $(SRC)/%.c
	mkdir -p '$(@D)'
	$(CC) -c $(INC) $(CFLAGS) $< -o $@






# todo make complete creating two targets
TARGETS = long_message_attack verify_hash



all: long_message_attack

# remove all $(TARGETS) members from dependencies
# and add long_message_attack.o as a dependency
long_message_attack: $(OBJDIR)/long_message_attack.o $(filter-out $(addsuffix .o, $(TARGETS)), $(OBJECTS) )
	$(CC)  $? $(LDFLAGS) $(LDLIBS) -o $@ 


.PHONY: clean
clean:
	rm -f $(OBJECTS)
	rm -f long_message_attack
