UNAME_S := $(shell uname -s)

CC = mpicxx
EXEC = exec
CCFLAGS = -std=c++11
INC = $(shell mpicc --showme:compile)
LDFLAGS = $(shell mpicc --showme:link)

SRC_FILES = $(wildcard src/*.cpp)
H_FILES = $(wildcard src/*.hpp)
OBJS = $(SRC_FILES:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(INC) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp %.hpp
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(EXEC) $(OBJS) *~ src/*~
