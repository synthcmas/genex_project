ifortErr = $(shell which ifort >/dev/null; echo $$?)
ifeq "$(ifortErr)" "0"
F90C     = ifort
FFLAGS   = -qopenmp -O3 -fpp -standard-semantics -no-wrap-margin 
DEBUGFLAGS= -qopenmp -O0 -ip -W0 -WB -fpp -traceback -check all -debug all -standard-semantics -no-wrap-margin
else
#gfortran compiler
F90C     = gfortran
FFLAGS   = -Ofast -fopenmp -ffree-form -ffree-line-length-none
DEBUGFLAGS = -cpp -Og -g -fbounds-check -fbacktrace -ffree-line-length-none -ffree-line-length-none
endif


# This is default code
Driver = driver
OBJFILES = rksuite.o parameters.o syseq.o solver.o utils.o inifile.o $(Driver).o

export FFLAGS
export DEBUGFLAGS
export F90C


%.o: %.f90
	$(F90C) $(FFLAGS) -c $*.f90
#	$(F90C) $(DEBUGFLAGS) -c $*.f90

all :	genex

genex:	$(OBJFILES)
	$(F90C) $(FFLAGS) $(OBJFILES) -o genex

debug:	$(OBJFILES)
	$(F90C) $(DEBUGFLAGS) $(OBJFILES) -o genex

run :
	time ./genex parameters.ini

plots :
	gnuplot > load 'plots.p'

clean:
	rm -f *.o *.mod *.d *.pc *.obj fort.* *.out plot_data/* genex 

