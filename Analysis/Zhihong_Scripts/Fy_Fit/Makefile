#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = FyFitMinuit
#PROGRAM = NkPlot
#PROGRAM = XSMinuitFit
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall
INCLUDE := -I${ROOTSYS}/include -I./SRC/ -I./
RTLIBS := -L${ROOTSYS}/lib -lCore -lMinuit -lgfortran 

all: $(PROGRAM)

GFORTRAN := gfortran
#GFORTRAN := f77
#FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

$./SRC/F1F209.o: ./SRC/F1F209.f 
	rm -f $@
	cd $./SRC/; $(GFORTRAN) $(FFLAGS) -c F1F209.f;


$(PROGRAM): $(PROGRAM).o ./SRC/F1F209.o
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) ./SRC/F1F209.o

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o ./SRC/*.o	
