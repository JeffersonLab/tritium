#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Check_XS
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)
XEMCDIR = /work/halla/e08014/disk1/yez/XEMC/XEMC_Clean/SRC

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include -I$(XEMCDIR)
RTLIBS := -L${ROOTSYS}/lib -lCore -lgfortran # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore ${x2}/mylib/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore /usr/lib64/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore  /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
all: $(PROGRAM)

GFORTRAN := gfortran
#GFORTRAN := f77
#FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

${XEMCDIR}/F1F209.o: ${XEMCDIR}/F1F209.f 
	rm -f $@
	cd ${XEMCDIR} ; $(GFORTRAN) $(FFLAGS) -c F1F209.f;

$(PROGRAM): $(PROGRAM).o $(XEMCDIR)/F1F209.o
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) $(XEMCDIR)/F1F209.o 

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o ./SRC/*.o	
