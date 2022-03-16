# Please run the following commands before make
# echo "export IRVSPDATA="`pwd` >> ~/.bashrc
# source ~/.bashrc

#ROOTDIR=../
#include ../Makefile
F90 = ifort
FLAGS = -g -traceback -check all -fpp 

REALMAKEFILE=Makefile.2 

default:lib

lib:
	$(MAKE) -f $(REALMAKEFILE) libs 

clean:
	rm -r *.mod *.o
