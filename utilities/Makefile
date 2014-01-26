#
# Makefile - very simple for now
#

#ifndef CXX
CXX = g++ -Wall -O3
#endif

CXXFLAGS = -g -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS

LIBOBJECTS = Table.o #RandomNumbers_TOS.o RandomNumbers_TNG.o
OBJECTS = $(LIBOBJECTS)

TARGETS = libUtility.a
LIBS = -lUtility

all: $(TARGETS)

static: LDFLAGS += =static
static: all

libUtility.a: $(LIBOBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) $(TARGETS) $(OBJECTS) *~

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

