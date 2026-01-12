#
# Makefile - for testing of the ebl models.
#

# Detect OS
UNAME_S := $(shell uname -s)

# Set compiler and flags based on OS
ifeq ($(UNAME_S),Darwin)
    # Note: Relies on variable `MACOS_DEPLOYMENT_TARGET` being set in the environment
    # export it as needed, e.g., `export MACOS_DEPLOYMENT_TARGET=14.7`
    CXX = clang++ -g -Wall -O3
    CXXFLAGS += $(INCFILES) -Xclang -fopenmp
    LIBS += $(EBLLIB) $(UTILITYLIB) -lomp
else ifeq ($(UNAME_S),Linux)
    CXX = g++ -g -Wall -O3
    CXXFLAGS += $(INCFILES) -fopenmp
    LIBS += $(EBLLIB) $(UTILITYLIB)
endif

# Directories
UTILITYDIR=../utilities
EBLDIR=../ebl

# Libraries
UTILITYLIB = -L$(UTILITYDIR) -lUtility
EBLLIB = -L$(EBLDIR) -lEBL

# Headers
INCFILES = -I$(EBLDIR) -I$(UTILITYDIR)

TARGETS = makeOptDepthTable

all: $(TARGETS)

makeOptDepthTable: makeOptDepthTable.o
	$(CXX) -o $@ $^ $(LIBS)

clean:
	$(RM) $(TARGETS) *~ *.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<