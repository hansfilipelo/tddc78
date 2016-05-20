UNAME=$(shell uname -s)

CC = mpic++
EXEC = physics
CCFLAGS = --std=c++14 -g -D_SIMULATION_STEPS_=100000
INC =
LDFLAGS = -lm

SRC_FILES = $(wildcard src/*.cpp)
H_FILES = $(wildcard src/*.hpp)
OBJS = $(SRC_FILES:.cpp=.o)

# Mac specific -----------

ifeq ($(UNAME), Darwin)
	INC += $(shell mpic++ --showme:compile)
	LDFLAGS += $(shell mpic++ --showme:link)
endif

# End Mac specific ------

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(INC) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp %.hpp
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(EXEC) $(OBJS)
