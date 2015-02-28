# Makefile for QSD code.

QSD_DIR = .

PYTHONPATH = /usr/local/usp/petens/python_install

CXX = g++
#INCLUDE = /opt/compiler/intel/ict/3.2/mkl/10.1/include
#INCLUDE = /usr/local/usp/COST/fftw-mpi-3.3.3/include
INCLUDE = $(PYTHONPATH)/include

#FLAGS = -g -O -Wno-deprecated  -I$(INCLUDE)
FLAGS = -g -O -Wno-deprecated  -Wno-write-strings -I$(INCLUDE)

#FC       = ifort  -c
FC       = gfortran  -c
#DFLAGS   =  -D__INTEL -D__LAPACK -D__FFTW3
DFLAGS   =  -D__LAPACK -D__FFTW3
##-D__FFTSG -D__FFTW2 -D__parallel -D__BLACS -D__SCALAPACK
#FCFLAGS =  -O2 -ftz -IPF-fp-relaxed  -heap-arrays 64 -I$(INCLUDE) 
FCFLAGS =  -O2 -I$(INCLUDE) 

INC = $(QSD_DIR)/include
SRC = $(QSD_DIR)/src

LOADLIBES = -lm
#LDFLAGS = -L/opt/compiler/intel/ict/3.2/mkl/10.1/lib/em64t -lmkl_scalapack

#LDFLAGS = -L$(MKLPATH) -lmkl_scalapack -lmkl_solver_ilp64_sequential -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -L/opt/compiler/intel/compiler/11.0/lib/intel64/ -lifcore $(LOADLIBES)

#LDFLAGS = -Wl,--start-group -Wl,--end-group -lpthread $(LOADLIBES)
LDFLAGS = -lgfortran -lpthread

LDFLAGS += $(PYTHONPATH)/lib/liblapack.a $(PYTHONPATH)/lib/libblas.a

#FFTFLAGS =  -L/usr/cta/unsupported/fftw/3.1.2-intel/lib -lfftw3
#FFTFLAGS =  -L/usr/local/usp/COST/fftw-mpi-3.3.3/lib -lfftw3
#FFTFLAGS =  /usr/local/usp/COST/fftw-mpi-3.3.3/lib/libfftw3.a
#FFTFLAGS =  -L/lustre/home1/u/richied/work2/lib -lfftw3
FFTFLAGS =  -L$(PYTHONPATH)/lib -lfftw3

#CXXFLAGS = $(DCC) -I$(INC) -I$(HOME)/include  $(FLAGS)
CXXFLAGS = $(DCC) -I$(INC) $(FLAGS)

#CPPFLAGS = -L$(MKLPATH) -lmkl_lapack -lmkl -lguide
CPPFLAGS = 

COMPILE = $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)
FCOMPILE  = $(FC) $(FCFLAGS) $(DFLAGS)
#LINK = icpc $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(FFTFLAGS)
#LINK = g++ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(FFTFLAGS)
LINK = g++ $(CPPFLAGS) $(CXXFLAGS) 
LOADLIBES += $(LDFLAGS) $(FFTFLAGS) -lm

RAN_OBJ =  Random.o Normal.o Uniform.o RNG.o ACG.o MLCG.o CmplxRan.o
RAN_HEAD = $(INC)/Random.h $(INC)/Normal.h $(INC)/Uniform.h $(INC)/RNG.h $(INC)/ACG.h $(INC)/MLCG.h $(INC)/CmplxRan.h

OBJ1 =  Operator.o PrimOp.o State.o Complex.o 
HEAD1 = $(INC)/Operator.h $(INC)/PrimOp.h $(INC)/State.h $(INC)/Complex.h $(INC)/Mesh.h
OBJ2 = FieldOp.o SpinOp.o AtomOp.o
OBJ3 = Traject.o $(RAN_OBJ)
OBJ4 = mesh.o coeff.o poisson.o 
OBJ5 = Mesh.o Coeff.o Epot.o Poisson.o
OBJ6 = cf.o fftw3.o
OBJ = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
#-------------------------------------------------------------------
# Add driver routines here.

ALL = onespin spins simple moving sums testprog template lineob damped qcascade sechar qd 

all: $(ALL)

qd: qd.cc $(OBJ)
	$(LINK) -o qd qd.cc $(OBJ)  $(LOADLIBES)

##qdf: qdf.f90 $(OBJ)
##	$(LINK) -o qdf qdf.f90 $(OBJ) $(LOADLIBES) 

onespin: onespin.cc $(OBJ)
	$(LINK) -o onespin onespin.cc $(OBJ) $(LOADLIBES)

spins: spins.cc $(OBJ)
	$(LINK) -o spins spins.cc $(OBJ) $(LOADLIBES)

simple: simple.cc $(OBJ)
	$(LINK) -o simple simple.cc $(OBJ) $(LOADLIBES)

lineob: lineob.cc $(OBJ)
	$(LINK) -o lineob lineob.cc $(OBJ) $(LOADLIBES)

qcascade: qcascade.cc $(OBJ)
	$(LINK) -o qcascade qcascade.cc $(OBJ) $(LOADLIBES)

sechar: sechar.cc $(OBJ)
	$(LINK) -o sechar sechar.cc $(OBJ) $(LOADLIBES)

damped: damped.cc $(OBJ)
	$(LINK) -o damped damped.cc $(OBJ) $(LOADLIBES)

moving: moving.cc $(OBJ)
	$(LINK) -o moving moving.cc $(OBJ) $(LOADLIBES)

sums: sums.cc $(OBJ)
	$(LINK) -o sums sums.cc $(OBJ) $(LOADLIBES)

testprog: testprog.cc $(OBJ)
	$(LINK) -o testprog testprog.cc $(OBJ) $(LOADLIBES)

template: template.cc $(OBJ)
	$(LINK) -o template template.cc $(OBJ) $(LOADLIBES)

#-------------------------------------------------------------------

cleanexe:
	-rm -f $(ALL)

cleanrand:
	-rm -f $(RAN_OBJ)

clean:
	-rm -f *.o

distclean: clean cleanrand cleanexe
	

Traject.o: $(SRC)/Traject.cc $(INC)/Traject.h $(HEAD1) $(RAN_HEAD)
	$(COMPILE) $(SRC)/Traject.cc

State.o: $(SRC)/State.cc $(HEAD1)
	$(COMPILE) $(SRC)/State.cc

Operator.o: $(SRC)/Operator.cc $(HEAD1)
	$(COMPILE) $(SRC)/Operator.cc

PrimOp.o: $(SRC)/PrimOp.cc $(HEAD1)
	$(COMPILE) $(SRC)/PrimOp.cc

FieldOp.o: $(SRC)/FieldOp.cc $(INC)/FieldOp.h $(HEAD1)
	$(COMPILE) $(SRC)/FieldOp.cc

SpinOp.o: $(SRC)/SpinOp.cc $(INC)/SpinOp.h $(HEAD1)
	$(COMPILE) $(SRC)/SpinOp.cc

AtomOp.o: $(SRC)/AtomOp.cc $(INC)/AtomOp.h $(HEAD1)
	$(COMPILE) $(SRC)/AtomOp.cc

Complex.o: $(SRC)/Complex.cc $(INC)/Complex.h
	$(COMPILE) $(SRC)/Complex.cc

Mesh.o: $(SRC)/Mesh.cc $(INC)/Mesh.h
	$(COMPILE) $(SRC)/Mesh.cc

Epot.o: $(SRC)/Epot.cc $(INC)/Epot.h
	$(COMPILE) $(SRC)/Epot.cc

Coeff.o: $(SRC)/Coeff.cc 
	$(COMPILE) $(SRC)/Coeff.cc

Poisson.o: $(SRC)/Poisson.cc
	$(COMPILE) $(SRC)/Poisson.cc

Random.o: $(SRC)/Random.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/Random.cc

Normal.o: $(SRC)/Normal.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/Normal.cc

Uniform.o: $(SRC)/Uniform.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/Uniform.cc

RNG.o: $(SRC)/RNG.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/RNG.cc

ACG.o: $(SRC)/ACG.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/ACG.cc

MLCG.o: $(SRC)/MLCG.cc $(RAN_HEAD)
	$(COMPILE) $(SRC)/MLCG.cc

CmplxRan.o: $(SRC)/CmplxRan.cc $(RAN_HEAD) $(INC)/Complex.h
	$(COMPILE) $(SRC)/CmplxRan.cc

coeff.o: $(SRC)/coeff.f90
	$(FCOMPILE) $(SRC)/coeff.f90

fftw3.o: $(SRC)/fftw3.f90
	$(FCOMPILE) $(SRC)/fftw3.f90

cf.o: $(SRC)/cf.f90
	$(FCOMPILE) $(SRC)/cf.f90

poisson.o: $(SRC)/poisson.f90
	$(FCOMPILE) $(SRC)/poisson.f90

mesh.o : $(SRC)/mesh.f90
	$(FCOMPILE) $(SRC)/mesh.f90
