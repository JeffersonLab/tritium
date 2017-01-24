#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Radiated
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include -I$(XEMDIR) -I$(ANALYZER)/src -I./
RTLIBS := -L${ROOTSYS}/lib -lCore -L${ANALYZER} -lHallA -ldc -lscaler /usr/lib/libg2c.so.0
all: $(PROGRAM)

$(PROGRAM): $(PROGRAM).o $(XEMDIR)/F1F209.o
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) $(XEMDIR)/F1F209.o 

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o	
