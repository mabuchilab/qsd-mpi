# Makefile for QSD code.
# Note: you may compile for a specific known platform with e.g.
#
#     make all PLATFORM=copper
#
#     or
#
#     make all PLATFORM=mlhpc
#

QSD_DIR = .
INC = $(QSD_DIR)/include
SRC = $(QSD_DIR)/src


# Compiler/linker flags common to C++ and Fortran
# -I : folder where C-compiler finds *.h Fortran compiler finds *.mod
# -O : optimization setting
# you could add debug options like -g here
FLAGS = -O2 -I$(QSD_DIR) -I$(INC)


# The C++ compiler
CXX        = g++
CXX_copper = CC
CXX_mlhpc  = mpiCC
ifdef PLATFORM
    CXX = ${CXX_$(PLATFORM)}
endif
COMPILE = $(CXX) -c $(FLAGS)


# The Fortran compiler
#FC       = ftn # copper
FC        = gfortran
FC_copper = ftn
FC_mlhpc  = gfortran
ifdef PLATFORM
    FC = ${FC_$(PLATFORM)}
endif
FCOMPILE  = $(FC) -c $(FLAGS)


# The Linker
# call to linker will be $(LINK) [objectfiles] $(LDFLAGS) $(LIBS)
# QSD depends on Lapack/BLAS and FFTW3
LINK = $(CXX) $(FLAGS) # Usually, the compiler executable knows how to link
LIBS        = -lblas -llapack -lgfortran -lm -lfftw3
LIBS_copper = -lfftw3 # note: `module load fftw`; LAPACK/BLAS are included in compiler
LIBS_mlhpc  = -lblas -llapack -lgfortran -lm -lmpi -lfftw3
ifdef PLATFORM
    LIBS = ${LIBS_$(PLATFORM)}
endif
LDFLAGS = -L$(QSD_DIR) # where .o files are located


#-------------------------------------------------------------------

# List all object files of the QSD library

RAN_OBJ =  Random.o Normal.o Uniform.o RNG.o ACG.o MLCG.o CmplxRan.o
RAN_HEAD = $(INC)/Random.h $(INC)/Normal.h $(INC)/Uniform.h $(INC)/RNG.h $(INC)/ACG.h $(INC)/MLCG.h $(INC)/CmplxRan.h

OBJ1 =  Operator.o PrimOp.o State.o Complex.o
HEAD1 = $(INC)/Operator.h $(INC)/PrimOp.h $(INC)/State.h $(INC)/Complex.h $(INC)/Mesh.h

OBJ2 = FieldOp.o SpinOp.o AtomOp.o
OBJ3 = Traject.o $(RAN_OBJ)
#OBJ4 = mesh.o coeff.o poisson.o
OBJ4 = mesh.o coeff.o
OBJ5 = Mesh.o Coeff.o Epot.o
OBJ6 = fftw3.o cube_function.o Poisson.o
OBJ = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)


#-------------------------------------------------------------------
# Add driver routines here.

ALL = onespin spins simple moving sums testprog template lineob damped qcascade sechar qd

all: $(ALL)

debug:
	$(warning *** SUMMARY OF VARIABLES ***)
	$(warning LINK = $(LINK))
	$(warning OBJ = $(OBJ))
	$(warning LDFLAGS = $(LDFLAGS))
	$(warning ****************************)

objectfiles: $(OBJ)

qd: qd.cc $(OBJ)
	$(LINK) -o qd qd.cc $(OBJ) $(LDFLAGS) $(LIBS)

onespin: onespin.cc $(OBJ)
	$(LINK) -o onespin onespin.cc $(OBJ) $(LDFLAGS) $(LIBS)

spins: spins.cc $(OBJ)
	$(LINK) -o spins spins.cc $(OBJ) $(LDFLAGS) $(LIBS)

simple: simple.cc $(OBJ)
	$(LINK) -o simple simple.cc $(OBJ) $(LDFLAGS) $(LIBS)

lineob: lineob.cc $(OBJ)
	$(LINK) -o lineob lineob.cc $(OBJ) $(LDFLAGS) $(LIBS)

qcascade: qcascade.cc $(OBJ)
	$(LINK) -o qcascade qcascade.cc $(OBJ) $(LDFLAGS) $(LIBS)

sechar: sechar.cc $(OBJ)
	$(LINK) -o sechar sechar.cc $(OBJ) $(LDFLAGS) $(LIBS)

damped: damped.cc $(OBJ)
	$(LINK) -o damped damped.cc $(OBJ) $(LDFLAGS) $(LIBS)

moving: moving.cc $(OBJ)
	$(LINK) -o moving moving.cc $(OBJ) $(LDFLAGS) $(LIBS)

sums: sums.cc $(OBJ)
	$(LINK) -o sums sums.cc $(OBJ) $(LDFLAGS) $(LIBS)

testprog: testprog.cc $(OBJ)
	$(LINK) -o testprog testprog.cc $(OBJ) $(LDFLAGS) $(LIBS)

template: template.cc $(OBJ)
	$(LINK) -o template template.cc $(OBJ) $(LDFLAGS) $(LIBS)

#-------------------------------------------------------------------

Traject.o: $(SRC)/Traject.cc $(INC)/Traject.h $(HEAD1) $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/Traject.cc

State.o: $(SRC)/State.cc $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/State.cc

Operator.o: $(SRC)/Operator.cc $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/Operator.cc

PrimOp.o: $(SRC)/PrimOp.cc $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/PrimOp.cc

FieldOp.o: $(SRC)/FieldOp.cc $(INC)/FieldOp.h $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/FieldOp.cc

SpinOp.o: $(SRC)/SpinOp.cc $(INC)/SpinOp.h $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/SpinOp.cc

AtomOp.o: $(SRC)/AtomOp.cc $(INC)/AtomOp.h $(HEAD1)
	$(COMPILE) -o $@ $(SRC)/AtomOp.cc

Complex.o: $(SRC)/Complex.cc $(INC)/Complex.h
	$(COMPILE) -o $@ $(SRC)/Complex.cc

Mesh.o: $(SRC)/Mesh.cc $(INC)/Mesh.h
	$(COMPILE) -o $@ $(SRC)/Mesh.cc

Epot.o: $(SRC)/Epot.cc $(INC)/Epot.h
	$(COMPILE) -o $@ $(SRC)/Epot.cc

Coeff.o: $(SRC)/Coeff.cc
	$(COMPILE) -o $@ $(SRC)/Coeff.cc

Poisson.o: $(SRC)/Poisson.cc
	$(COMPILE) -o $@ $(SRC)/Poisson.cc

Random.o: $(SRC)/Random.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/Random.cc

Normal.o: $(SRC)/Normal.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/Normal.cc

Uniform.o: $(SRC)/Uniform.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/Uniform.cc

RNG.o: $(SRC)/RNG.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/RNG.cc

ACG.o: $(SRC)/ACG.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/ACG.cc

MLCG.o: $(SRC)/MLCG.cc $(RAN_HEAD)
	$(COMPILE) -o $@ $(SRC)/MLCG.cc

CmplxRan.o: $(SRC)/CmplxRan.cc $(RAN_HEAD) $(INC)/Complex.h
	$(COMPILE) -o $@ $(SRC)/CmplxRan.cc

coeff.o: $(SRC)/coeff.f90
	$(FCOMPILE) -o $@ $(SRC)/coeff.f90

fftw3.o: $(SRC)/fftw3.f90
	$(FCOMPILE) -o $@ $(SRC)/fftw3.f90

cube_function.o: $(SRC)/cube_function.f90
	$(FCOMPILE) -o $@ $(SRC)/cube_function.f90

poisson.o: $(SRC)/poisson.f90
	$(FCOMPILE) -o $@ $(SRC)/poisson.f90

mesh.o: $(SRC)/mesh.f90
	$(FCOMPILE) -o $@ $(SRC)/mesh.f90


cleanexe:
	rm -f $(ALL)

cleanrand:
	rm -f $(RAN_OBJ)

clean:
	rm -f *.o *.mod

distclean: clean cleanexe

#-------------------------------------------------------------------

.PHONY: all debug objectfiles cleanexe cleanrand clean distclean

