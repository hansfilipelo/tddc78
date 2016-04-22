
CC = mpicc
BLUR_EXEC = blur
THRES_EXEC = thres
CCFLAGS = -std=c99 -DMAX_RAD=1000 -DMAX_PIXELS=1000*1000
INC = $(shell mpicc --showme:compile)
LDFLAGS = $(shell mpicc --showme:link)

BLUR_SRC_FILES = $(wildcard src/blur*.c)
BLUR_SRC_FILES += src/ppm.c src/ppmio.c src/gaussw.c src/timing_mach.c
BLUR_H_FILES = $(wildcard src/blur*.h) $(wildcard src/ppm*.h)
BLUR_H_FILES += src/ppm.h src/ppmio.h src/gaussw.h
BLUR_OBJS = $(BLUR_SRC_FILES:.c=.o)

$(BLUR_EXEC): $(BLUR_OBJS)
	$(CC) $(BLUR_OBJS) $(INC) -o $(BLUR_EXEC) $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CCFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(BLUR_EXEC) $(BLUR_OBJS) 
