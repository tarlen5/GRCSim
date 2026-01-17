/*!

  makeOptDepthTable.cpp

  Makes a table (to use with Table2D class) of:

    -1    z1               z2                z3...
    egy1 tau(egy1,z1) tau(egy1,z2)....
    egy2 ...
    ...

  and saves it to an output file of type:
    optDepth_EBLModel<>_z<>_0.01TeV_1TeV.txt

*/

#include "DIRBR.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void InitializeEBL(string &ebl_modelname, DIRBR &ebl) {

  ifstream model;
  model.open(ebl_modelname.c_str());

  // if (model == NULL) {
  if (!model.is_open()) {
    std::cerr << "ERROR: could not open EBL model file: " << ebl_modelname
              << std::endl;
    exit(EXIT_FAILURE);
  }

  DIRBR_ERR ebl_err;
  vector<double> lambda_vec;
  vector<double> nuFnu_vec;
  // Set up DIRBR curve...
  while (model) {
    double lambda_i;
    double nuFnu_i;
    model >> lambda_i >> nuFnu_i;
    if (model) {
      lambda_vec.push_back(lambda_i);
      nuFnu_vec.push_back(nuFnu_i);
    }
  }
  model.close();

  // Initialize DIRBR:
  ebl_err = ebl.SetDIRBR(lambda_vec.size(), &lambda_vec[0], &nuFnu_vec[0]);
  std::cout << "SetDIRBR returned: " << ebl_err << std::endl;
  ebl_err = ebl.TestDIRBRlimits();
  std::cout << "TestDIRBRlimits returned: " << ebl_err << std::endl;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void DefineBinning(vector<double> &egyBinsVec, vector<double> &zBinsVec,
                   double egyMin, double egyMax, double redshift) {

  // double egyMin = 0.01;   // TeV
  // double egyMax = 10.0;  // TeV
  double egyBinsDec = 4.0;
  double egy = egyMin;
  double ibin = 0.0;
  while (egy <= egyMax) {
    egy = egyMin * pow(10.0, ibin / egyBinsDec);
    egyBinsVec.push_back(egy);
    ibin += 1.0;
  }

  double zMin = 0.01;
  double zMax = redshift;
  double zBinsDec = 4.0;
  double z = zMin;
  ibin = 0.0;
  while (z <= zMax) {
    z = zMin * pow(10.0, ibin / zBinsDec);
    zBinsVec.push_back(z);
    ibin += 1.0;
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int main(int argc, char **argv) {

  char *progname = *argv;
  argc--, argv++;

  // Check for help option
  if (argc == 1 &&
      (std::string(argv[0]) == "-h" || std::string(argv[0]) == "--help")) {
    std::cout
        << "Usage: " << progname
        << " EBLModel redshift egy_low egy_high [num_threads]\n\n"
        << "Arguments:\n"
        << "  EBLModel     : Name of the EBL model file (e.g., EBLModel4msld)\n"
        << "  redshift     : Maximum redshift value (e.g., 0.3)\n"
        << "  egy_low      : Minimum energy in TeV (e.g., 0.01)\n"
        << "  egy_high     : Maximum energy in TeV (e.g., 10.0)\n"
        << "  num_threads  : Optional number of threads for parallel "
           "computation (default: 1)\n\n"
        << "Example: " << progname << " EBLModel4msld 0.3 0.01 10.0 4\n"
        << "This generates optDepth_EBLModel4msld_z0.3_0.01TeV_10.0TeV.txt\n";
    exit(0);
  }

  // Simple get options:
  if (argc < 4 || argc > 5) {
    std::cerr << "ERROR! Requires 4 or 5 arguments. Use -h or --help for usage "
                 "information.\n"
              << "Example: " << progname
              << " EBLModel4msld 0.3 0.01 10.0 [num_threads]\n";
    exit(EXIT_FAILURE);
  }
  string eblmodel = *argv;
  argv++;

  string s_z = *argv;
  istringstream ss_z(s_z);
  argv++;
  double redshift;
  ss_z >> redshift;

  string s_egyLo = *argv;
  istringstream ss_elo(s_egyLo);
  argv++;
  double egyLo = 0.0;
  ss_elo >> egyLo;

  string s_egyHi = *argv;
  istringstream ss_ehi(s_egyHi);
  argv++;
  double egyHi = 0.0;
  ss_ehi >> egyHi;

  int num_threads = 1;
  if (argc == 5) {
    string s_threads = *argv;
    istringstream ss_threads(s_threads);
    ss_threads >> num_threads;
  }
  omp_set_num_threads(num_threads);
  ///////////////////////////////////////////////////////

  string outfilename = "optDepth_" + eblmodel + "_z" + s_z + "_" + s_egyLo +
                       "TeV_" + s_egyHi + "TeV_" + to_string(num_threads) +
                       "cores.txt";

  /////////////////////////////////////////////////////////////
  /////////// Set up Egy, z binning   /////////////////////////
  /////////////////////////////////////////////////////////////
  vector<double> egyBinsVec;
  vector<double> zBinsVec;
  DefineBinning(egyBinsVec, zBinsVec, egyLo, egyHi, redshift);

  string ebl_modelfile = "../" + eblmodel + ".dat";
  DIRBR ebl;
  InitializeEBL(ebl_modelfile, ebl);

  ofstream outfile(outfilename.c_str());
  outfile << "-1 ";

  for (unsigned iz = 0; iz < zBinsVec.size(); iz++)
    outfile << zBinsVec[iz] << " ";
  outfile << endl;

  std::cout << "Calculating optical depth table..." << endl;
  std::cout << "[tau], [energy TeV], [redshift]" << endl;
  for (unsigned iegy = 0; iegy < egyBinsVec.size(); iegy++) {
    double egy = egyBinsVec[iegy];
    outfile << egy << " ";
    vector<double> taus(zBinsVec.size());

#pragma omp parallel for
    for (unsigned iz = 0; iz < zBinsVec.size(); iz++) {
      double redshift = zBinsVec[iz];
      taus[iz] = ebl.OpticalDepth(egy, redshift);
#pragma omp critical
      std::cout << taus[iz] << " " << egy << " " << redshift << std::endl;
    }
    for (auto tau : taus)
      outfile << tau << " ";
    outfile << endl;
  }
  outfile.close();

  return 0;
}
