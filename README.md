GRCSim
======

Monte Carlo Simulation of gamma-ray initiated electromagnetic cascades in intergalactic space.

Relies on two externals:
1. the qd library, a high-precision software library developed at LBNL, and can be obtained here: http://crd-legacy.lbl.gov/~dhbailey/mpdist/
2. the cern ROOT library, which can be obtained here: http://root.cern.ch/drupal/content/downloading-root

Once these are up and running, some adjustments to the makefiles must be made, since I haven't been able to automate this process yet. Please email me at: tca3@psu.edu for details and permission to use and modify this code.

TODO:

1. Make high precision propagation its own class/module so code common to GalacticGrid and PairProduction classes are using same function from include file.
2. Implement a sqlite database solution for magnetic field files, rather than storing them as text files.
