#  makefile

# -----------------------------------------------------------------------------
# Root
# -----------------------------------------------------------------------------
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`

#--------------------------------------
# lib and include directories
#--------------------------------------
# qd/dd directory:
includeqddir=qd_dd/qd-2.3.13/include
INCLUDE=include/
UTILDIR = utilities/
EBLDIR   = ebl/

# LIB DIRS
#LIBQDDIR = qd_dd/lib/
LIBLOCAL=lib/

# OBJECT DIR:
OBJDIR=objdir/

#SRC DIR:
SRC=source/

# Compiler
CC=gcc
CPP=g++
CFLAGS= -g -Wall -O3 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -Wgnu-static-float-init
#Note when debugging, use -O0 optimization flag...

#-----------Use this version if running gprof-----------
#CPP=g++ -pg
#CFLAGS= -Wall -O0 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
#-------------------------------------------------------

CFLAGS += -I$(includeqddir) -I$(INCLUDE) -I$(UTILDIR) -I$(EBLDIR)

# Archive
AR=ar

# common libraries
LIBS= $(ROOTLIBS) -lProp -lqd -lpthread -lEBL -lUtility
LDFLAGS= -L$(LIBLOCAL) -L$(EBLDIR) -L$(UTILDIR)

TARGETS = $(OBJDIR) run_sim_cascade
#OBJ = anyoption.o IGCascadeSim.o PairProduction.o KleinNishina.o \
#	  RelParticle.o Vec4D.o Vec3D.o GalacticGrid.o
OBJ = $(subst $(SRC), $(OBJDIR), $(patsubst %.cpp, %.o, $(wildcard $(SRC)*.cpp)))

all: $(TARGETS)

run_sim_cascade: run_sim_cascade.o $(LIBLOCAL)libProp.a
	$(CPP) $(LDFLAGS) -o $@ $< $(LIBS)
	mv $< $(OBJDIR)

$(LIBLOCAL)libProp.a: $(OBJ)
	$(AR) -rc $(LIBLOCAL)libProp.a $(OBJ)

$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	$(RM) -r lib/*.a $(TARGETS) *~ *.o

$(OBJDIR)%.o: $(SRC)%.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $< -o $@

%.o: %.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $<

#%.o: source/%.cpp
#	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $<

#%.o: %.cpp
#	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $<
