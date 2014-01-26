#
# Makefile
#


#ifndef CXX
CXX           = g++ -Wall -O3
#endif


UTILITYDIR = ../utilities/

# Libraries
UTILITY = -L$(UTILITYDIR) -lUtility

#CXXFLAGS += 
INCFILES = -I$(UTILITYDIR)

# Add them all up
CXXFLAGS = -g $(INCFILES) -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
LIBS = $(UTILITY)

LIBOBJECTS = DIRBR.o DIRBR_Redshift.o
OBJECTS = $(LIBOBJECTS)

TARGETS = libEBL.a

LIBS = -lEBL

all: $(TARGETS)

static: LDFLAGS += -static
static: all

libEBL.a: $(LIBOBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) $(TARGETS) $(OBJECTS) *~

%.o: $.cpp
	$(CXX) $(CXXFLAGS) -c $<
