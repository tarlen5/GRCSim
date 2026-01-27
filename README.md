GRCSim
======

Monte Carlo simulation of gamma-ray initiated electromagnetic cascades in intergalactic space. Simulation code written in collaboration with others in the gamma ray astrophysics group at UCLA from 2004 - 2013. (Revised for ease of use since then).

<h2> Installation Notes </h2>

Relies on two external dependencies:

1. The qd library, which can be obtained from [this repository](https://github.com/BL-highprecision/QD). To install the qd library, read the instructions to set up autoconfig on your system, then perform the following steps:
    ```
    ./configure
    make
    make install
    ```
  to install into the standard places. After installation, to link during compliation and execution of GRCSim, define:
  - `QDLIB=/usr/local/lib` (or wherever the file `libqd.so` is located after install)
  - `QDINC=/usr/local/include` (or wherever the file `qd/dd_real.h` is located after install)

2. The HDF5 library for storing data, which can be installed on most operating systems easily.


Finally, modify your LD_LIBRARY_PATH as such:
  - `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QDLIB`

Once these are installed, and environment variables set, if you are on a Mac/Linux system, the only thing needed to compile is to run `make`. Please email me at: timothyarlen@gmail.com for questions, further information, and permission to use and modify this code.

<h2> Running the executable to perform a Cascade Simulation</h2>

An executable is created called `run_sim_cascade`. In order to view the help and all options, type:

  ```./run_sim_cascade -h```

For example, to run one iteration of a $1$ TeV gamma photon cascade at a redshift of $0.14$ with a Intergalactic Magnetic Field strength of $10^{-12}$ G and coherence length of $1$ Mpc, run:

  `./run_sim_cascade 1000 1e-12 1 0.14 0 1`

Before the simulation is able to run, it needs access to two input files:
  1. Optical Depth Files for the ebl model used. This can be created from `OptDepthFiles/makeOptDepthTable.cpp`. See instructions there for how to create the optical depth files at the redshift and gamma photon energies desired.
  2. Magnetic Field File - If running for the first time, create the magnetic field file directory and filename needed (Or the code will create the directory and touch an empty file for you - NOT ADVISED when running on computing cluster due to possibile race conditions).
