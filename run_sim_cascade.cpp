/*!
---------------------------------------------------------------------------
  \file     run_sim_cascade
            run_sim_cascade.cpp

      Produces a single gamma photon cascade in intergalactic space,
      a number of times, as determined by the user.

      This file is basically a wrapper to call the IGCascadeSim class,
      which performs the intergalactic cascade simulation, with the
      given parameters.

  \author   Timothy C. Arlen
            timothyarlen@gmail.com

  \date     Jan 24, 2014

  \note     Initial code I'm putting under version control, after major
            rewrite of input variables, options, configurations, etc. into
      a much more streamlined and easier version using the AnyOption
      class.
---------------------------------------------------------------------------
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "IGCascadeSim.hpp"
#include "anyoption.h"

using namespace std;

AnyOption *DefineOptions(int argc, char *argv[], const string &progname) {

  AnyOption *opt = new AnyOption();

  /* 1. SET THE USAGE/HELP   */
  opt->addUsage("Usage:\n");
  string usageHelp = "  " + progname +
                     " energy Bmag Lcoh ze num_gammas file_num" +
                     " [Options]\n";
  opt->addUsage(usageHelp.c_str());
  opt->addUsage("energy     - [GeV] energy of photon at OBSERVER redshift.");
  opt->addUsage("Bmag       - [gauss] magnetic field magnitude.");
  opt->addUsage("Lcoh       - [Mpc] coherence length of magnetic field.");
  opt->addUsage("ze         - redshift at source point.");
  opt->addUsage("file_num   - file number tag.");
  opt->addUsage("iterations - number of photons to simulate in this run.\n");
  // opt->addUsage( "Usage: " );
  opt->addUsage("Options:");
  opt->addUsage(" -h  --help      Prints this help ");
  opt->addUsage(" --eblmodel      EBLModel4msld        EBLModel name.");
  opt->addUsage(
      " --mf_dir        MagneticFieldFiles/  magnetic field files dir."
  );
  opt->addUsage(" --opt_depth_dir OptDepthFiles/       opt depth files dir.");
  opt->addUsage(" --output_dir    (cwd)      Directory to output files.");
  opt->addUsage(" --gam_egy_min    0.1 min energy [GeV] to track gammas.");
  opt->addUsage(" --lep_egy_min    75. min energy [GeV] to track leptons.");
  opt->addUsage(
      " --seed           <uint32_t> Random number seed for reproducibility."
  );

  opt->addUsage(" --mf_no_lock    No locking of mf files, so cosmic variance "
                "is not preserved in simulation. ");
  opt->addUsage(" --single_gen    Forces single generation of cascade.");
  opt->addUsage(" --trk_delay     Track time delay throughout cascade.");
  opt->addUsage(" --trk_leptons   Tracks all leptons throughout cascade.");

  /* 2. SET THE OPTION STRINGS/CHARACTERS */
  opt->setFlag("help", 'h');
  opt->setOption("eblmodel");
  opt->setOption("mf_dir");
  opt->setOption("opt_depth_dir");
  opt->setOption("output_dir");
  opt->setOption("gam_egy_min");
  opt->setOption("lep_egy_min");
  opt->setOption("seed");
  opt->setFlag("mf_no_lock");
  opt->setFlag("single_gen");
  opt->setFlag("trk_delay");
  opt->setFlag("trk_leptons");

  return opt;
}

int main(int argc, char **argv) {

  string progname(*argv);
  AnyOption *opt = DefineOptions(argc, argv, progname);

  // Process All
  opt->processCommandArgs(argc, argv);
  if (!opt->hasOptions()) {
    opt->printUsage();
    exit(EXIT_FAILURE);
  }

  // Access values:
  if (opt->getFlag("help") || opt->getFlag('h')) {
    opt->printUsage();
    exit(0);
  }

  if (opt->getArgc() != 6) {
    cerr << "\nERROR: 6 mandatory inputs required! Usage: " << endl;
    opt->printUsage();
    exit(EXIT_FAILURE);
  }

  // Process Command Line Args:
  string egy = opt->getArgv(0);
  string mag_field = opt->getArgv(1);
  string coh_len = opt->getArgv(2);
  string redshift = opt->getArgv(3);
  string file_count = opt->getArgv(4);
  string iterations = opt->getArgv(5);

  cout << "\n>> Arguments parsed: " << endl;
  cout << "  egy:        " << egy << endl;
  cout << "  redshift:   " << redshift << endl;
  cout << "  mag_field:  " << mag_field << endl;
  cout << "  coh_len:    " << coh_len << endl;
  cout << "  file_count: " << file_count << endl;
  cout << "  iterations: " << iterations << endl;

  IGCascade::IGCascadeSim *my_sim = new IGCascade::IGCascadeSim(
      egy, redshift, mag_field, coh_len, file_count, opt
  );

  int num_iterations = atoi(iterations.c_str());
  my_sim->RunCascade(num_iterations);

  return 0;
}
