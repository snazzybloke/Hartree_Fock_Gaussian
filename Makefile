FC=ifort
FFLAGS= -O3

EXEC = abmol

MKLROOT=/opt/intel/oneapi/mkl/latest
MKLPATH=$(MKLROOT)/lib/intel64

INCLUDE=$(MKLROOT)/include
INC=-I$(INCLUDE) -I$(INCLUDE)/intel64/lp64

BLAS=-Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_lapack95_lp64.a -Wl,--end-group

LAPACK=-L$(MKLPATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64
NAG=-lpthread -L/home/ante/lib/ -lnag

# statically linked:
# LIBS= $(BLAS) $(NAG) 
# dynamically linked:
LIBS= $(LAPACK) $(NAG)

OBJS = abtype.o koordtype.o basis_sets.o store_integrals.o matrix.o abmol.o integrals.o collector.o scf.o ext_utils.o

%.o: %.f90
	$(FC) -c $(FFLAGS)  $<  $(INC)


abmol: $(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LIBS) 

clean:
	rm *.mod *.o
