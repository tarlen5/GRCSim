#  makefile

SHELL=/bin/sh
MAKE=/usr/bin/make
ECHO=echo

# -----------------------------------------------------------------------------
# Root
# -----------------------------------------------------------------------------
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`

#--------------------------------------
# lib and include directories
#--------------------------------------
# qd/dd directory:
#QDLIB = qd_dd/build_hammer_2014_Nov_13_2.3.13/lib
QDLIB = /usr/local/lib
UTILDIR = utilities/
EBLDIR   = ebl/

# Archive
AR=ar

#LOCAL DIRS:
INCLUDE=include/
LIBLOCAL=lib/
OBJDIR=objdir/
SRC=source/

# Compiler
CC=gcc
CPP=g++
CFLAGS= -g -Wall -O3
# -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS

#-----------Use this version if running gprof/gdb-----------
#CPP=g++ -pg
#CFLAGS= -Wall -O0 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
#-----------------------------------------------------------

CFLAGS += -I$(QDINC) -I$(INCLUDE) -I$(UTILDIR) -I$(EBLDIR)

# common libraries
LIBS= $(ROOTLIBS) -lProp -lqd -lpthread -lEBL -lUtility
LDFLAGS= -L$(LIBLOCAL) -L$(EBLDIR) -L$(UTILDIR) -L$(QDLIB)
OBJLIBS=libUtility.a libEBL.a

TARGETS = $(OBJDIR) $(OBJLIBS) run_sim_cascade
OBJ = $(subst $(SRC), $(OBJDIR), $(patsubst %.cpp, %.o, $(wildcard $(SRC)*.cpp)))

.PHONY: $(OBJLIBS)

all: $(TARGETS)

# run_sim_cascade: run_sim_cascade.o $(LIBLOCAL)libProp.a $(OBJLIBS)
# 	$(CPP) $(LDFLAGS) -o $@ $< $(LIBS)
# 	mv $< $(OBJDIR)

run_sim_cascade: run_sim_cascade.o $(LIBLOCAL)libProp.a $(OBJLIBS)
	$(CPP) $(LDFLAGS) $< $(LIBS) $(ROOTLIBS) -Iinclude/ -o $@
	mv $< $(OBJDIR)

$(LIBLOCAL)libProp.a: $(OBJ)
	$(AR) -rc $(LIBLOCAL)libProp.a $(OBJ)

libUtility.a:
	$(MAKE) -C $(UTILDIR)

libEBL.a:
	$(MAKE) -C $(EBLDIR)

$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	$(RM) -r lib/*.a $(TARGETS) *~ *.o
	$(RM) -rf $(OBJDIR)
	$(MAKE) clean -C $(UTILDIR)
	$(MAKE) clean -C $(EBLDIR)

$(OBJDIR)%.o: $(SRC)%.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $< -o $@ -Iinclude/

%.o: %.cpp
	$(CPP) $(CFLAGS) $(ROOTCFLAGS) -c $< -Iinclude/

