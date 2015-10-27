GRCSim
======

Monte Carlo Simulation of gamma-ray initiated electromagnetic cascades in intergalactic space.

Relies on two externals:

1. The qd library, a high-precision software library developed at LBNL, and can be obtained here: http://crd-legacy.lbl.gov/~dhbailey/mpdist/
	1. May need an older compiler to compile this software (worked for me on gcc v. 4.1.2)
	2. You will need to note the paths to the library and include directories for compliling the GRCSim code (see below)
2. The cern ROOT library, which can be obtained here: http://root.cern.ch/drupal/content/downloading-root

ENVIRONMENT VARIABLES:
* In order to compile the GRCSim code, the following environment variables must be set: 
  * Path to the installed qd library: QDLIB=/path/to/qdlibdir
  * Path to qd library include files: QDINC=/path/to/qdinclude

Once these are installed, and environment variables set, if you are on a Mac/Linux system, the only thing needed to compile is to run `make`. Please email me at: tca3@psu.edu for questions, further information, and permission to use and modify this code.

TODO:

1. Implement an sqlite database solution for magnetic field files, rather than storing them as text files, to account for cosmic variance over a specific simulation run.
