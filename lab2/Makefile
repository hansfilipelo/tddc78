HOSTNAME=$(shell hostname)

CF=gfortran
CFFLAGS= -fopenmp #-DOMP_NUM_THREADS=4
EXEC=laplsolv

ifeq ($(HOSTNAME), triolith1)
	CF=ifort
endif


all: $(EXEC)

$(EXEC): src/laplsolv.f90
	$(CF) $(CFFLAGS) $^ -o $@

clean:
	rm -rf $(EXEC)
