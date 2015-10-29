#
# Makefile - for testing of the ebl models.
#

CXX = g++ -g -Wall -O3

# Directories
UTILITYDIR=../utilities
EBLDIR=../ebl

# Libraries
UTILITYLIB = -L$(UTILITYDIR) -lUtility
EBLLIB = -L$(EBLDIR) -lEBL

# Add them all up
#CXXFLAGS += 

# Headers
INCFILES = -I$(EBLDIR) -I$(UTILITYDIR)

#Add up
CXXFLAGS += $(INCFILES)
LIBS += $(EBLLIB) $(UTILITYLIB)

TARGETS = makeOptDepthTable

all: $(TARGETS)

makeOptDepthTable: makeOptDepthTable.o
	$(CXX) -o $@ $^ $(LIBS)

clean:
	$(RM) $(TARGETS) *~ *.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
