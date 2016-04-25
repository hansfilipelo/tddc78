
CC = mpicc
BLUR_EXEC = blur
THRES_EXEC = thres

CCFLAGS = -std=c99 -DMAX_RAD=10000 -DMAX_PIXELS=MAX_RAD*MAX_RAD
INC = $(shell mpicc --showme:compile)
LDFLAGS = $(shell mpicc --showme:link)
COMMON_FILES = src/ppm.c src/ppmio.c src/timing_mach.c src/mpi_data_types.c
COMMON_H_FILES = src/ppm.h src/ppmio.h src/mpi_data_types.h

BLUR_SRC_FILES = $(wildcard src/blur*.c) src/gaussw.c $(COMMON_FILES)
BLUR_H_FILES = $(wildcard src/blur*.h) $(wildcard src/ppm*.h) src/gaussw.h $(COMMON_H_FILES)
BLUR_OBJS = $(BLUR_SRC_FILES:.c=.o)

THRES_SRC_FILES = $(wildcard src/thres*.c) $(COMMON_FILES)
THRES_H_FILES = $(wildcard src/thres*.h) $(wildcard src/ppm*.h) $(COMMON_H_FILES)
THRES_OBJS = $(THRES_SRC_FILES:.c=.o)

all: $(THRES_EXEC) $(BLUR_EXEC)

$(THRES_EXEC): $(THRES_OBJS)
	$(CC) $(THRES_OBJS) $(INC) -o $(THRES_EXEC) $(LDFLAGS)

$(BLUR_EXEC): $(BLUR_OBJS)
	$(CC) $(BLUR_OBJS) $(INC) -o $(BLUR_EXEC) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(BLUR_EXEC) $(BLUR_OBJS) $(THRES_EXEC) $(THRES_OBJS)
