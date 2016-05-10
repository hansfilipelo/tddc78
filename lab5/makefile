UNAME=$(shell uname -s)

CC = mpicc
EXEC = physics
CCFLAGS = --std=gnu99 -DMAX_RAD=10000 -DMAX_PIXELS=MAX_RAD*MAX_RAD
INC =
LDFLAGS = -lm

SRC_FILES = $(wildcard src/*.c)
H_FILES = $(wildcard src/*.h)
OBJS = $(SRC_FILES:.c=.o)

# Mac specific -----------

ifeq ($(UNAME), Darwin)
	INC += $(shell mpicc --showme:compile)
	LDFLAGS += $(shell mpicc --showme:link)
	COMMON_FILES += src/timing_mach.c
	COMMON_H_FILES += src/timing_mach.h
endif

# End Mac specific ------

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(INC) -o $(EXEC) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(EXEC) $(OBJS)
