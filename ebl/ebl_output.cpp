/*!

  make_ebl_output.cpp - Tests the new DIRBR class: DIRBR_Redshift.cpp,
    which includes a redshift evolution ebl model (in this case, Dominguez
    2012), and plots lambda I_lambda vs lambda.

 */

#include "DIRBR.hpp"
#include "DIRBR_Redshift.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void InitializeEBL(DIRBR *dirbr, string eblmodelfile);

int main(int argc, char **argv) {

  char *progname = *argv;
  argc--, argv++;

  if (argc != 1) {
    cerr << "ERROR! Please use 1 input. Ex: \n  " << progname
         << " EBLModel.dat\n";
  }

  string eblmodelfile = *argv;
  argv++;

  if (true) {
    DIRBR_Redshift *dirbr = new DIRBR_Redshift(eblmodelfile);

    string output = "output_redshift_" + eblmodelfile;
    ofstream outfile(output.c_str());
    double lambda0 = 0.1; // [mum]
    double bins_per_dec = 16.0;
    double lambda = 0.0;
    double ibin = 0.0;
    double z = 1.0;
    while (lambda < 1000.0) {
      lambda = lambda0 * pow(10.0, ibin / bins_per_dec);
      double value = dirbr->GetDIRBR(lambda, z);
      outfile << lambda << " " << value << endl;
      ibin++;
    }
    outfile.close();

  } else {
    DIRBR *dirbr = new DIRBR();
    InitializeEBL(dirbr, eblmodelfile);

    string output = "output_redshift_" + eblmodelfile;
    ofstream outfile(output.c_str());
    double lambda0 = 0.1; // [mum]
    double bins_per_dec = 16.0;
    double lambda = 0.0;
    double ibin = 0.0;
    while (lambda < 1000.0) {
      lambda = lambda0 * pow(10.0, ibin / bins_per_dec);
      double value = dirbr->GetDIRBR(lambda, 0.0);
      outfile << lambda << " " << value << endl;
      ibin++;
    }
    outfile.close();
  }

  return 0;
}

void InitializeEBL(DIRBR *dirbr, string eblmodelfile) {

  std::ifstream model;
  model.open(eblmodelfile.c_str());

  if (model == NULL) {
    std::cerr << "ERROR: could not open EBL model file: " << eblmodelfile
              << std::endl;
    exit(EXIT_FAILURE);
  }

  vector<double> m_lambda_vec;
  vector<double> m_nuFnu_vec;
  DIRBR_ERR m_ebl_err;
  // Set up DIRBR curve...
  while (model) {
    double lambda_i;
    double nuFnu_i;
    model >> lambda_i >> nuFnu_i;
    if (model) {
      m_lambda_vec.push_back(lambda_i);
      m_nuFnu_vec.push_back(nuFnu_i);
    }
  }
  model.close();

  // Initialize DIRBR:
  m_ebl_err =
      dirbr->SetDIRBR(m_lambda_vec.size(), &m_lambda_vec[0], &m_nuFnu_vec[0]);
  std::cout << "SetDIRBR returned: " << m_ebl_err << std::endl;
  m_ebl_err = dirbr->TestDIRBRlimits();
  std::cout << "TestDIRBRlimits returned: " << m_ebl_err << std::endl;
}
