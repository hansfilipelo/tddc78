HOSTNAME=$(shell hostname)

CF=gfortran
CFFLAGS= -DOMP_NUM_THREADS=4
EXEC=laplsolv

ifeq ($(HOSTNAME), triolith1)
	CF=ifort
	CCFLAGS += -fopenmp
else
	CCFLAGS += -fopenmp
endif


all: $(EXEC)

$(EXEC): src/laplsolv.f90
	$(CF) $(CCFLAGS) $^ -o $@

clean:
	rm -rf $(EXEC)

