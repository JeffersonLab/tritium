#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

  PROGRAM = VDC_EffL
# PROGRAM = VDC_EffR
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall
INCLUDE := -I${ROOTSYS}/include -I$(XEMDIR)
RTLIBS := -L${ROOTSYS}/lib -lCore -L${ANALYZER} -lHallA -ldc -lscaler

all: $(PROGRAM)

$(PROGRAM): $(PROGRAM).o
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS)

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o	
