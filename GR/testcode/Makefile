.IGNORE =
.SUFFIXES: .F90 .f90

all: bssn_headers eqtest trans


####### intel compilers
#F90_COMMAND    = ifort
#F90FLAGS       = -O3 -r8
#F90FIXEDFORM   = -FI
#F90FREEFORM    = -FR

#######  GNU compilers
CXX_COMMAND    = g++-6 -std=c++14
CXXFLAGS       = -g -march=native -ftree-vectorize -fgcse-after-reload -finline-functions -funswitch-loops
F90_COMMAND    = gfortran
#F90FLAGS       = -g -fdefault-real-8 -ffpe-trap=invalid,zero
F90FLAGS       = -g -fdefault-real-8 
F90FIXEDFORM   = -ffixed-form -ffixed-line-length-132
F90FREEFORM    = -ffree-form

INCFLAGS       = -I./include -I/usr/local/include
#LIBS           = -L$(HOME)/gnu/lib -lbbhutil -lrt -lgfortran
LIBS           = -L/usr/local/lib -lgfortran

F90            = $(F90_COMMAND) $(F90FLAGS) $(INCFLAGS)
CXX            = $(CXX_COMMAND) $(CXXFLAGS) $(INCFLAGS)
CXX_LOAD       = $(CXX_COMMAND)

OBJ   = main.o bssnrhs.o rhs.o rhsopt.o optderivs.o hadrhs.o hadbssnrhs.o rhsvec.o

.F90.o:
	$(F90) $(F90FIXEDFORM) -c $*.F90

.f90.o:
	$(F90) $(F90FREEFORM) -c $*.f90

.cc.o:
	$(CXX) -c $*.cc


bssn_headers:
	cd include && $(MAKE)

eqtest: $(OBJ)
	$(CXX_LOAD) $(OBJ) $(LIBS) -o eqtest

trans: opt_gradients.o
	$(CXX_LOAD) $^ -o $@
	
check_adv: check_adv.o
	$(CXX) $^ -o $@


clean:
	cd include && $(MAKE) clean
	/bin/rm -f *.o *.mod *.optrpt eqtest trans 
	/bin/rm -f *sdf

