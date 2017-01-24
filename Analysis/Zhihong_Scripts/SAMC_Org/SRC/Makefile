# hamc = Hall A Monte Carlo
# R. Michaels, May 2008
# 
# Modified for compiling SAMC
# Z. Ye, August 2011

# Choose the compiler.
GCC=g++
GLD=g++

MAKENODICTIONARY=1

export OSNAME := $(shell uname)
export ROOTSYS_32=/u/apps/root/5.26-00/root
export ROOTSYS_64=/u/apps/root/5.26-00-64bit/root

ifeq ($(OSNAME),SunOS)

   ROOTCFLAGS    = $(shell root-config --cflags)
   ROOTLIBS      = $(shell root-config --libs)
   ROOTGLIBS     = $(shell root-config --glibs)
   CXX           = $(GCC)
   CXXFLAGS      = -KPIC -DSUNVERS -I$(ROOTSYS_32)/include -I$(MAINDIR) #-I$(XEMDIR)
#   CXXFLAGS      = -KPIC -DSUNVERS -I$(ROOTSYS_64)/include -I$(MAINDIR) -I$(XEMDIR)
   CXXFLAGS     += $(ROOTCFLAGS)
   LD            = $(GLD)
   LDFLAGS       = -g -D
   SOFLAGS       = -G
   GLIB =  -lm -lc -lgen -lw -lnsl -ldl
   SLIB = -L/opt/SUNWspro/SC4.2/lib -lGFORTRAN -lM77 -lsunmath
   ET_AC_FLAGS = -D_REENTRANT -D_POSIX_THREAD_SEMANTICS
   ET_CFLAGS = -mt -fast -xO5 -KPIC $(ET_AC_FLAGS) -DSUNVERS
   ONLIBS = -lposix4 -lnsl -lsocket -lresolv
   LIBS = $(GLIB)

endif

# Linux with egcs

ifeq ($(OSNAME),Linux)
   ROOTLIBS      = $(shell root-config --libs)
   ROOTGLIBS     = $(shell root-config --glibs)
   INCLUDES      = -I$(ROOTSYS_32)/include -I$(XEMDIR)
 #  INCLUDES      = -I$(ROOTSYS_64)/include -I$(XEMDIR)
   CXX           = $(GCC)
   CXXFLAGS      = -Wall -fno-exceptions -fPIC $(INCLUDES)
   LD            = $(GLD)
   LDFLAGS       = 
   SOFLAGS       = -shared 
   GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 
   ET_AC_FLAGS = -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
   ET_CFLAGS = -02 -fPIC -I. $(ET_AC_FLAGS) -DLINUXVERS
   ONLIBS = -lieee -lpthread -ldl -lresolv 
   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) -lgfortran -lgfortranbegin -lnsl
   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) /usr/lib/libg2c.so.0 
#   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS)  /home/yez/mylib/libg2c.so.1
endif

MAKEDEPEND    = $(GCC)

ALL_LIBS = $(LIBS) 

SRCDIR=./
INCLUDES += -I$(SRCDIR)

ifdef PROFILE
   CXXFLAGS += -pg
endif

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g
endif

  CXXFLAGS += -g	

SRC = SAMC.C 
#withxycollimator.C


DEPS = $(SRC:.C=.d)
DEP  = $(SRC:.C=.d)
HEAD = $(SRC:.C=.h) 

PROGS = SAMC

# Make the dictionary
ifdef MAKENODICTIONARY
  OBJS = $(SRC:.C=_NODICT.o)
else
  OBJS = $(SRC:.C=.o)
endif

##### ######

GFORTRAN := f77
#FFLAGS = -c -g -Wall -Wextra -Wconversion -pedantic
#FFLAGS = -ffixed-line-ffixed-line-length-0 -x f77-cpp-input-length-0 -x f77-cpp-input  \
         -fno-automatic -ff2c -fsecond-underscore -fbounds-check
FFLAGS = -C -g  -fbounds-check -finit-local-zero -fno-automatic 

GFORTRAN_DIR =  ./
GFORTRAN_SRC =  $(GFORTRAN_DIR)/monte_trans_hrs.f $(GFORTRAN_DIR)/Left_funcs.f $(GFORTRAN_DIR)/Right_funcs.f \
           $(GFORTRAN_DIR)/Left_r-function.f $(GFORTRAN_DIR)/Right_r-function.f
GFORTRAN_OBJ =   $(GFORTRAN_DIR)/monte_trans_hrs.o $(GFORTRAN_DIR)/Left_funcs.o $(GFORTRAN_DIR)/Right_funcs.o \
           $(GFORTRAN_DIR)/Left_r-function.o $(GFORTRAN_DIR)/Right_r-function.o $(XEMDIR)/F1F209.o
GFORTRAN_OBJS += $(GFORTRAN_OBJ) 


install: all 

all: $(PROGS) ${GFORTRAN_OBJS}

SAMC:  $(OBJS) $(SRC) $(HEAD) $(GFORTRAN_OBJS)
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(GFORTRAN_OBJS) $(ALL_LIBS)

withxycollimator.o: withxycollimator.C
	rm -f $@
	$(LD) $(CXXFLAGS) -c $@  withxycollimator.C $(OBJS) $(GFORTRAN_OBJS) $(ALL_LIBS);

$(GFORTRAN_DIR)/monte_trans_hrs.o: $(GFORTRAN_DIR)/monte_trans_hrs.f 
	rm -f $@
	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c monte_trans_hrs.f;
$(GFORTRAN_DIR)/Left_funcs.o: $(GFORTRAN_DIR)/Left_funcs.f 
	rm -f $@
	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Left_funcs.f;
$(GFORTRAN_DIR)/Right_funcs.o: $(GFORTRAN_DIR)/Right_funcs.f 
	rm -f $@
	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Right_funcs.f;

$(GFORTRAN_DIR)/Right_r-funcs.o: $(GFORTRAN_DIR)/Right_r-funcs.f 
	rm -f $@
#	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Right_r-funcs.f;
	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Right_r-funcs.f;
$(GFORTRAN_DIR)/Left_r-funcs.o: $(GFORTRAN_DIR)/Left_r-funcs.f 
	rm -f $@
#	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Left_r-funcs.f;
	cd $(GFORTRAN_DIR) ; $(GFORTRAN) $(FFLAGS) -c Left_r-funcs.f;


$(XEMDIR)/F1F209.o: $(XEMDIR)/F1F209.f 
	rm -f $@
	cd $(XEMDIR) ; $(GFORTRAN) $(FFLAGS) -c F1F209.f;


#############################   XGT2 OBJS   ###############

tarfile: clean version
	tar cvf $(VERS).tar ./$(VERS)

clean:
	rm -f $(GFORTRAN_DIR)/*.o 

realclean:  clean
	rm -f *.d *.tar  *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%_NODICT.o:	%.C
	$(CXX) $(CXXFLAGS) -c -DNODICT -o $*_NODICT.o $<

%.d:	%.C 
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
                | sed '\''s%^.*\.o%$*\.o%g'\'' \
                | sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
                [ -s $@ ] || rm -f $@'

-include $(DEPS)
