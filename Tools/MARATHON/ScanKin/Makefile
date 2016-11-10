#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

#PROGRAM = ScanBBB
PROGRAM = BB_Xb2
#PROGRAM = ScanBB
#PROGRAM = BB_Xb
#PROGRAM = BB_Ep
#PROGRAM = BB_Plot
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include
RTLIBS := -L${ROOTSYS}/lib -lCore #-lHistPainter 
all: $(PROGRAM)

$(PROGRAM): $(PROGRAM).o 
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS)

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o 
