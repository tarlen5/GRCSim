#  makefile

# -----------------------------------------------------------------------------
# Root
# -----------------------------------------------------------------------------
ROOTCONFIG=root-config
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTGLIBS:= $(shell $(ROOTCONFIG) --glibs)

#--------------------------------------
# lib and include directories
#--------------------------------------
# qd/dd directory:
includeqddir=qd_dd/qd-2.3.13/include
includeusr=include/
UTILDIR = utilities/
EBLDIR   = ebl/

# LIB DIRS
#LIBQDDIR = qd_dd/lib/
LIBLOCAL = lib/

# Compiler
CC=gcc
CPP=g++
CFLAGS= -g -Wall -O3 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -Wgnu-static-float-init
#Note when debugging, use -O0 optimization flag...

#-----------Use this version if running gprof-----------
#CPP=g++ -pg	
#CFLAGS= -Wall -O0 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
#-------------------------------------------------------

CFLAGS += -I$(includeqddir) -Iinclude/ -I$(UTILDIR) -I$(EBLDIR)

# Archive
AR=ar

# common libraries
LIBS= $(ROOTGLIBS) -lProp -lqd -lpthread -lEBL -lUtility
#LDFLAGS= -L$(LIBLOCAL) -L$(LIBQDDIR) -L$(EBLDIR) -L$(UTILDIR)
LDFLAGS= -L$(LIBLOCAL) -L$(EBLDIR) -L$(UTILDIR)

TARGETS = run_sim_cascade
OBJECTS = anyoption.o IGCascadeSim.o PairProduction.o KleinNishina.o \
	  RelParticle.o Vec4D.o Vec3D.o GalacticGrid.o

all: $(TARGETS)

run_sim_cascade: run_sim_cascade.o $(LIBLOCAL)libProp.a
	$(CPP) $(LDFLAGS) -o $@ $< $(LIBS)

$(LIBLOCAL)libProp.a: $(OBJECTS) 
	$(AR) -rc $(LIBLOCAL)libProp.a $(OBJECTS)

clean:
	$(RM) lib/*.a $(TARGETS) *~ *.o

%.o: source/%.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $<

%.o: %.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $<
