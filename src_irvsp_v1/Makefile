.SUFFIXES: .o .f .f90
F90    = ifort 
LINKER = $(F90) 
FLAGS =   -g -traceback -check all #-heap-arrays 2048 

MODULES=  init.o symm.o 
WAVEMD =  wave_data.o
VECT = irrep.o 
VECT2= wrtdate.o latgen2.o rmprop.o raxis.o su2op.o eulang.o invmat.o invmati.o \
       kgroup.o pntgrp.o dprodci.o crwcnd.o classe.o trsyma.o rotkv.o chrct.o wrtir.o \
       mdfpg.o  trsymb.o
ALL = $(MODULES) $(WAVEMD) $(VECT) $(VECT2)


default:irvsp
	
irvsp:$(ALL)
	$(LINKER) $(FLAGS) -o irvsp $(ALL)  $(LIBS)
	#rm -f *.o *.mod

$(MODULES):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(WAVEMD):%.o:%.f90
	ifort -assume byterecl -c  $*.f90
$(VECT):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(VECT2):%.o:%.f
	$(F90) -c -FR  $*.f

clean:
	rm -f *.mod *.o 
