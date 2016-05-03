HOSTNAME=$(shell hostname)

CF=gfortran
EXEC=laplsolv

ifeq ($(HOSTNAME), triolith1)
	CF=ifort
endif


all: $(EXEC)

$(EXEC): src/laplsolv.f90
	$(CF) -fopenmp $^ -o $@

clean:
	rm -rf $(EXEC)
