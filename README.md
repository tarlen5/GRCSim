GRCSim
======

Monte Carlo simulation of gamma-ray initiated electromagnetic cascades in intergalactic space. Simulation code written in collaboration with others in the gamma ray astrophysics group at UCLA from 2004 - 2013.

<h2> Installation Notes </h2>

Relies on two external dependencies:

1. The qd library, which can be obtained here: http://crd-legacy.lbl.gov/~dhbailey/mpdist/. To install the qd library, download the `.tar.gz` file from the above webpage and perform the following steps:
    ```
    ./configure
    make
    make install
    ```
  to install into the standard places. After installation, to link during compliation and execution of GRCSim, define:
  - `QDLIB=/usr/local/lib` (or wherever the file `libqd.so` is located after install)
  - `QDINC=/usr/local/include` (or wherever the file `qd/dd_real.h` is located after install)

2. The cern ROOT library, which can be obtained here: http://root.cern.ch/drupal/content/downloading-root. As of 18 Oct 2023, installing the most recent version of ROOT (6.28.06).


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
  2. Magnetic Field File - If running for the first time, create the magnetic field file directory and filename needed (instructions for that will be given in the error message when the executable, `run_sim_cascade` is run).

<h2>TODO</h2>

* Implement the Intergalactic Magnetic Field in a more realistsic way.
  1. Currently it uses a "GalacticGrid" class which models the Magnetic Field as constant within a 3D cell with edge length=coherence length of the magnetic field. The values are stored within a `MagneticFieldFiles/` directory which is a text file of the cell coordinates and magnetic field direction. Possibly an sqlite database solution for the magnetic field files would be better when running on a culster, rather than storing them as text files.
  2. Alternatively, maybe the whole 3D grid needs to be re-thought with another implementation.


