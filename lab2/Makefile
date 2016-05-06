UNAME=$(shell uname -s)

CC = gcc
EXEC = blur
CCFLAGS = --std=gnu99 -Wall -DMAX_RAD=1000 -DMAX_PIXELS=MAX_RAD*MAX_RAD -D_N_TASKS_=4
INC =
LDFLAGS = -lpthread

COMMON_FILES= src/mpi_data_types.c src/ppm.c src/ppmio.c
COMMON_H_FILES= src/ppm.h src/ppmio.h

# Mac specific -----------

ifeq ($(UNAME), Darwin)
	INC +=
	LDFLAGS +=
	COMMON_FILES += src/timing_mach.c src/barrier.c
	COMMON_H_FILES += src/timing_mach.h src/barrier.h
endif

# End Mac specific ------

SRC_FILES = $(wildcard src/blur*.c) $(COMMON_FILES) src/gaussw.c
H_FILES = $(wildcard src/blur*.h) $(COMMON_H_FILES)
OBJS = $(SRC_FILES:.c=.o)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(INC) -o $(EXEC) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(EXEC) $(OBJS) $(THRES_EXEC) $(THRES_OBJS)
