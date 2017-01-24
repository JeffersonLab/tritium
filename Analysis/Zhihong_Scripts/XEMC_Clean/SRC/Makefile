#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = Check_XS 
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall
INCLUDE := -I${ROOTSYS}/include -I${XEMDIR}
RTLIBS := -L${ROOTSYS}/lib -lCore /usr/lib/libg2c.so.0 
all: $(PROGRAM)

GFORTRAN := f77
FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 

GFORTRAN_OBJ =  ${XEMDIR}/F1F209.o
GFORTRAN_OBJS += $(GFORTRAN_OBJ) 

$(PROGRAM): $(PROGRAM).o $(GFORTRAN_OBJS)
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS) $(GFORTRAN_OBJS)

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

${XEMDIR}/F1F209.o: ${XEMDIR}/F1F209.f 
	rm -f $@
	cd ${XEMDIR} ; $(GFORTRAN) $(FFLAGS) -c F1F209.f;

clean:
#	rm -f *.o ${GFORTRAN_OBJ}
	rm -f *.o
