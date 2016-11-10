#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = RateCal
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)
XEMCDIR = /Users/yez/work/github/XEMC/SRC

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include -I$(XEMCDIR)
RTLIBS := -L${ROOTSYS}/lib -lCore -lgfortran 
all: $(PROGRAM)

GFORTRAN := gfortran
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

${XEMCDIR}/F1F209.o: ${XEMCDIR}/F1F209.f 
	rm -f $@
	cd ${XEMCDIR} ; $(GFORTRAN) $(FFLAGS) -c F1F209.f;

$(PROGRAM): $(PROGRAM).o $(XEMCDIR)/F1F209.o
	clang++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) $(XEMCDIR)/F1F209.o 

$(PROGRAM).o: $(PROGRAM).C
	clang++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o ./SRC/*.o ${XEMCDIR}/*.o	
