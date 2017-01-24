 #-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Ep_Yield
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

INCLUDE := -I${ROOTSYS}/include -I${XGT2_XS}
RTLIBS := -L${ROOTSYS}/lib -lCore -L${ANALYZER} -lHallA -ldc -lscaler
CFLAGS += -Wall -fno-exceptions -fPIC $(INCLUDES)
all: $(PROGRAM)

$(PROGRAM): $(PROGRAM).o
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS)
$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C

clean:
	rm -f *.o	
