.SUFFIXES: .o .f .f90
F90    = ifort  
LINKER = $(F90) 
FLAGS  = -g -traceback #-check all 

MODULES=  init.o symm.o
MOD1   =  nonsymm.o
WAVEMD =  wave_data.o
VECT = irrep.o addsign.o invreal33.o kgroup.o chrct.o 
VECT2= wrtdate.o latgen2.o rmprop.o raxis.o su2op.o eulang.o invmat.o invmati.o \
       pntgrp.o pntgrpb.o dprodci.o crwcnd.o classe.o trsyma.o rotkv.o \
       wrtir.o trsymb.o
ALL =  $(MOD1) $(MODULES) $(WAVEMD) $(VECT)  $(VECT2)

default: irvsp
	
tar:
	tar -cvzf tb.tar.gz *.f90  plot.gnu makefile *.data *.in
	
irvsp:$(ALL)
	$(LINKER) $(FLAGS) -o irvsp  $(ALL)  $(LIBS)
	
#https://www.ruanyifeng.com/blog/2015/02/make.html
#3.7 自动变量（Automatic Variables）

#Make命令还提供一些自动变量，它们的值与当前规则有关。主要有以下几个。

#（1）$@

#$@指代当前目标，就是Make命令当前构建的那个目标。比如，make foo的 $@ 就指代foo。


#    a.txt b.txt: 
#        touch $@

#等同于下面的写法。


#    a.txt:
#        touch a.txt
#    b.txt:
#        touch b.txt

#（2）$<

#$< 指代第一个前置条件。比如，规则为 t: p1 p2，那么$< 就指代p1。


#    a.txt: b.txt c.txt
#        cp $< $@ 

#等同于下面的写法。


#    a.txt: b.txt c.txt
#        cp b.txt a.txt 

#（3）$?

#$? 指代比目标更新的所有前置条件，之间以空格分隔。比如，规则为 t: p1 p2，其中 p2 的时间戳比 t 新，$?就指代p2。

#（4）$^

#$^ 指代所有前置条件，之间以空格分隔。比如，规则为 t: p1 p2，那么 $^ 就指代 p1 p2 。

#（5）$*

#$* 指代匹配符 % 匹配的部分， 比如% 匹配 f1.txt 中的f1 ，$* 就表示 f1。

#（6）$(@D) 和 $(@F)

#$(@D) 和 $(@F) 分别指向 $@ 的目录名和文件名。比如，$@是 src/input.c，那么$(@D) 的值为 src ，$(@F) 的值为 input.c。

#（7）$(<D) 和 $(<F)

#$(<D) 和 $(<F) 分别指向 $< 的目录名和文件名。

#所有的自动变量清单，请看手册。下面是自动变量的一个例子。


#    dest/%.txt: src/%.txt
#        @[ -d dest ] || mkdir dest
#        cp $< $@

#上面代码将 src 目录下的 txt 文件，拷贝到 dest 目录下。首先判断 dest 目录是否存在，如果不存在就新建，然后，$< 指代前置文件（src/%.txt）， $@ 指代目标文件（dest/%.txt）。

$(MOD1):%.o:%.f90
	$(F90) -c $(FLAGS) -fpp -DIRVSPDATA $*.f90

#$(MOD1):part1 part2
#	cat part1 >  nonsymm.f90
#	echo "      spgpath = '$(CURDIR)'" >> nonsymm.f90
#	cat part2 >> nonsymm.f90
#	$(F90) -c $(FLAGS) $*.f90
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
