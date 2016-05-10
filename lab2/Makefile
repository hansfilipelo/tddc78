UNAME=$(shell uname -s)

CC = gcc
BLUR_EXEC = blur
THRES_EXEC = thres
CCFLAGS = -g --std=gnu99 -Wall -DMAX_RAD=10000 -DMAX_PIXELS=MAX_RAD*MAX_RAD -D_N_TASKS_=4
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

ifeq ($(UNAME), Linux)
	LDFLAGS += -lm -lrt
endif

# End Mac specific ------

BLUR_SRC_FILES = $(wildcard src/blur*.c) $(COMMON_FILES) src/gaussw.c
BLUR_H_FILES = $(wildcard src/blur*.h) $(COMMON_H_FILES)
BLUR_OBJS = $(BLUR_SRC_FILES:.c=.o)

THRES_SRC_FILES = $(wildcard src/thres*.c) $(COMMON_FILES)
THRES_H_FILES = $(wildcard src/thres*.c) $(COMMON_H_FILES)
THRES_OBJS = $(THRES_SRC_FILES:.c=.o)

all: $(BLUR_EXEC) $(THRES_EXEC)

$(BLUR_EXEC): $(BLUR_OBJS)
	$(CC) $(BLUR_OBJS) $(INC) -o $(BLUR_EXEC) $(LDFLAGS)

$(THRES_EXEC): $(THRES_OBJS)
	$(CC) $(THRES_OBJS) $(INC) -o $(THRES_EXEC) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(BLUR_EXEC) $(BLUR_OBJS) $(THRES_EXEC) $(THRES_OBJS)
