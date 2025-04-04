#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
# FC       = ifort -assume byterecl 
FC       = mpiifort -assume byterecl 

FFLAGS   = -g -O3 -traceback  -check bounds  -qopenmp # -parallel -heap-arrays
MAKE     = make

MKL_PATH = $(MKLROOT)/lib/intel64
MKL_INC  = $(MKLROOT)/include
MKLFLAG  = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread
BLAS     = $(MKL_PATH)/libmkl_blas95_lp64.a
LAPACK   = $(MKL_PATH)/libmkl_lapack95_lp64.a
LLIBS    = $(BLAS) $(LAPACK) 

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC = prec.f90 utils.f90 parallel.f90 fileio.f90 couplings.f90 hamil.f90 \
   	TimeProp.f90 dish.f90 fssh.f90 main.f90


OBJ = $(SRC:.f90=.o)
EXE = hfnamd

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<
.SUFFIXES: .o .F90
.F90.o:
	$(FC) $(FFLAGS) -c $< -DENABLEMPI -DENABLEMKL

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
tdm:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) -I$(MKL_INC) -L$(MKL_PATH) $(MKLFLAG) 
# $(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(LLIBS)

clean:
	rm -f *.mod
	rm -f *.o
tar:
	tar -czvf hfnamd.tgz *.f90 Makefile
tag:
	ctags *.f90
