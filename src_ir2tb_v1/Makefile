.SUFFIXES: .o .f .f90
F90    = ifort   #-Wl, -stack_size,0x10000000 # -warn all
LINKER = $(F90) 
FLAGS =   -g -traceback -check all #-heap-arrays 2048 
LIBS = -L${MKLROOT}/lib  -mkl

MODULES=  init.o symm.o  wave_data.o
VECT = irrep.o  file_util.o
VECT2= wrtdate.o latgen2.o rmprop.o raxis.o su2op.o eulang.o invmat.o invmati.o \
       kgroup.o pntgrp.o dprodci.o crwcnd.o classe.o trsyma.o rotkv.o chrct.o wrtir.o \
       mdfpg.o  trsymb.o 
ALL = $(MODULES)  $(VECT)  $(VECT2)


default:irre
	
tar:
	tar -cvzf tb.tar.gz *.f90  plot.gnu makefile *.data *.in
	
irre:$(ALL)
	$(LINKER) $(FLAGS) -o ir2tb $(ALL)  $(LIBS)
	rm -f wave_data.mod

$(MODULES):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(VECT):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(VECT2):%.o:%.f
	$(F90) -c -FR  $*.f


run:
	mpirun -np 4  ./banding.exe
clean:
	rm -f *.mod *.o DTOT.data 
