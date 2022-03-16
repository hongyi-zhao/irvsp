.SUFFIXES: .o .f .f90
F90    = ifort 
LINKER = $(F90) 
FLAGS  =  -g -traceback #-check all #-heap-arrays 2048 
#LIBS   = -L${MKLROOT}/lib/em64t -mkl

MODULES=  init.o symm.o 
WAVEMD =  wave_data.o
VECT   =  irrep.o inv_int33.o kgroup.o rot_KGvec.o chrct.o wrtir.o
ALL    = $(MODULES) $(WAVEMD) $(VECT) 

default:irvsp
	
irvsp:$(ALL)
	$(LINKER) $(FLAGS) -o vasp2trace $(ALL)  $(LIBS)
	rm -f *.o *.mod

$(MODULES):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(WAVEMD):%.o:%.f90
	ifort -assume byterecl -c  $*.f90
$(VECT):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90

clean:
	rm -f *.mod *.o *.data 
