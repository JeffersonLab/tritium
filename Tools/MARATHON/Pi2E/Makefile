#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Pi2E

CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall
all: $(PROGRAM)

$(PROGRAM): $(PROGRAM).o 
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS)  $(LIBS) $(GLIBS)

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C 

clean:
	rm -f *.o
