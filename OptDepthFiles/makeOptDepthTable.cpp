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
#include <sstream>
#include <string>
#include <vector>
#include <fstream>


using namespace std;

void InitializeEBL(string& ebl_modelname, DIRBR& ebl)
{
  
  ifstream model;
  model.open(ebl_modelname.c_str());
  
  if (model == NULL) {
    std::cerr << "ERROR: could not open EBL model file: " << ebl_modelname 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  DIRBR_ERR ebl_err;
  vector<double> lambda_vec;
  vector<double> nuFnu_vec;
  // Set up DIRBR curve...
  while(model) {
    double lambda_i;
    double nuFnu_i;
    model >> lambda_i >> nuFnu_i;
    if(model) {
      lambda_vec.push_back(lambda_i);
      nuFnu_vec.push_back(nuFnu_i);
    }
  }
  model.close();
  
  // Initialize DIRBR:
  ebl_err=ebl.SetDIRBR(lambda_vec.size(), &lambda_vec[0], 
		       &nuFnu_vec[0]);
  std::cout << "SetDIRBR returned: " << ebl_err << std::endl;
  ebl_err=ebl.TestDIRBRlimits();
  std::cout << "TestDIRBRlimits returned: " << ebl_err << std::endl;
 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


void DefineBinning(vector<double>& egyBinsVec, vector<double>& zBinsVec,
		   double egyMin, double egyMax, double redshift)
{
  
  //double egyMin = 0.01;   // TeV
  //double egyMax = 10.0;  // TeV
  double egyBinsDec = 4.0;
  double egy = egyMin;
  double ibin = 0.0;
  while(egy <= egyMax) {
    egy = egyMin*pow(10.0,ibin/egyBinsDec);
    egyBinsVec.push_back(egy);
    ibin+=1.0;
  }
  
  double zMin   = 0.01;
  double zMax   = redshift;
  double zBinsDec = 4.0;
  double z = zMin;
  ibin = 0.0;
  while(z <= zMax) {
    z = zMin*pow(10.0,ibin/zBinsDec);
    zBinsVec.push_back(z);
    ibin+=1.0;
  }
  
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


int main(int argc, char** argv) {
  
  char* progname = *argv;
  argc--,argv++;

  // Simple get options:
  if( argc != 4) {
  cerr<<"ERROR! Please use 4 inputs. Ex: \n  "<<progname<<
    " EBLModel4msld redshift egy_low [TeV] egy_high [TeV] \n";
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
  ///////////////////////////////////////////////////////
    

  string outfilename = 
    "optDepth_"+eblmodel+"_z"+s_z+"_"+s_egyLo+"TeV_"+s_egyHi+"TeV.txt";
  

  /////////////////////////////////////////////////////////////
  /////////// Set up Egy, z binning   /////////////////////////
  /////////////////////////////////////////////////////////////
  vector<double> egyBinsVec;
  vector<double> zBinsVec;
  DefineBinning(egyBinsVec, zBinsVec, egyLo, egyHi, redshift);

  
  string ebl_modelfile = "/gpfs/home/tca3/work/GRCSim/"+eblmodel+".dat";
  DIRBR ebl;
  InitializeEBL(ebl_modelfile, ebl);

  ofstream outfile(outfilename.c_str());
  outfile<<"-1 ";
  
  for(unsigned iz = 0; iz<zBinsVec.size(); iz++)
    outfile<<zBinsVec[iz]<<" ";
  outfile<<endl;

  for(unsigned iegy = 0; iegy<egyBinsVec.size(); iegy++) {
    double egy = egyBinsVec[iegy];
    outfile<<egy<<" ";
    for(unsigned iz = 0; iz<zBinsVec.size(); iz++) {
      double redshift = zBinsVec[iz];

      double tau = ebl.OpticalDepth(egy,redshift);      
      //outfile << tau <<" "<<egy<<" "<<redshift<<endl;
      outfile<<tau<<" ";
      cout<< tau <<" "<<egy<<" "<<redshift<<endl;
    }
    outfile<<endl;
  }
  outfile.close();
  
  return 0;

}

