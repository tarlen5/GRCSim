/*!
-------------------------------------------------------------------------------
    \file   IGCascadeSim.hpp

    Header file for IGCascadeSim.cpp

    \author   Timothy C. Arlen                      \n
              Department of Physics and Astronomy   \n
              UCLA                                  \n
              arlen@astro.ucla.edu                  \n

    \date     May 20, 2012                          \n

    \note
-------------------------------------------------------------------------------
*/

#ifndef IGCASCADE_IGCASCADESIM_H
#define IGCASCADE_IGCASCADESIM_H

#include "DIRBR.hpp"
#include "GalacticGrid.hpp"
#include "HighPrecProp.hpp"
#include "KleinNishina.hpp"
#include "PairProduction.hpp"
#include "PhysicsConstants.hpp"
#include "RandomNumbers.hpp"
#include "RelParticle.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "anyoption.h"

// HDF5 INCLUDES:
#include <H5Cpp.h>
#include <hdf5.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <string>
#include <time.h>
#include <vector>

using namespace H5;
using namespace std;

namespace IGCascade {

class IGCascadeSim {
public:
  IGCascadeSim(
      const string &egy, const string &redshift, const string &mag_field,
      const string &coh_len, const string &file_count, AnyOption *opt
  );

  inline VEC3D_T GetMinLepEgy() { return m_egy_lepton_min; }

  void CreateDirectGamma(RelParticle &GammaPhoton);

  void RunCascade(const int numIterations);

private:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Functions   ///////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void ProcessOptions(
      AnyOption *opt, string &eblmodel, string &mf_dir, string &opt_depth_dir,
      string &output_dir
  );
  string DefineCascadeFile(
      const string &s_eblmodel, const string &s_egy, const string &s_Bmag,
      const string &s_ze, const string &s_cellsize, const string &s_file_num,
      const string &output_dir
  );
  string DefineMFfile(
      const string &static_var_file, const string &s_Bmag,
      const string &s_cellsize, const string &s_ze
  );
  string DefineLowEgyFile(
      const string &s_eblmodel, const string &s_egy, const string &s_Bmag,
      const string &s_ze, const string &s_cellsize, const string &s_file_num
  );
  string DefineTrackLeptonFile(
      const string &s_eblmodel, const string &s_egy, const string &s_Bmag,
      const string &s_ze, const string &s_cellsize, const string &s_file_num,
      const string &output_dir
  );
  string DefineTrackTimeDelayFile(
      const string &s_eblmodel, const string &s_egy, const string &s_Bmag,
      const string &s_ze, const string &s_cellsize, const string &s_file_num
  );
  void DefineRp(void);
  void InitializeEBL(const string &ebl_model_file);
  string DefineOptDepthTable(
      const string &opt_depth_dir, const string &s_eblmodel, const string &s_ze
  );
  void RunSinglePhotonCascade(void);
  void PropagateDirectPhoton(
      RelParticle &GammaPhoton, std::stack<RelParticle *> &lepton_stack
  );
  void DefineLeptons(
      RelParticle &GammaPhoton, RelParticle &EBLPhoton, RelParticle *Electron,
      RelParticle *Positron
  );
  void CreateLeptonTree(RelParticle *Lepton);
  void WriteLeptonToFile(RelParticle *Lepton);
  double GetOpticalDepthVal(double egy, double z);
  void PropagatePhotonToObserver(RelParticle &GammaPhoton);
  bool PropagateSecondaryPhoton(
      RelParticle &GammaPhoton, RelParticle *Electron, RelParticle *Positron
  );
  double GetMinZ(
      const Vec4D &gam_ph_p4, const Vec4D &gam_ph_r4, const VEC3D_T &gam_ph_z
  );
  void PropagateLepton(
      std::stack<RelParticle *> &lepton_stack, RelParticle &GammaPhoton
  );
  void StoreSecPhoton(Vec4D &gam_ph_p4, Vec4D &gam_ph_r4, const double &weight);
  void SaveToLowEnergyFile(RelParticle &Particle);
  void SaveDirectPhoton(RelParticle &Particle, std::ofstream &photon_list);
  void SaveLepton(RelParticle *Lepton);
  void SaveToTrackTimeDelayFile(RelParticle &Particle);

  // HDF5 writing functions:
  void WriteSecPhotonsToHDF5();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Data Members///////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Option variables:
  bool m_single_gen_bool;
  bool m_trk_delay_bool;
  bool m_trk_leptons_bool;
  VEC3D_T m_egy_gamma_min;
  VEC3D_T m_egy_lepton_min;
  bool m_LOCK;
  bool m_seed_provided;
  uint32_t m_seed_value;

  // Cascade parameters
  VEC3D_T m_egy_cascade; // TeV
  VEC3D_T m_bmag;
  VEC3D_T m_ze;
  string m_cellsize;

  unsigned m_globalLeptonNum;

  string m_cascade_file;
  string m_low_egy_file;
  string m_save_lepton_file;
  string m_track_lepton_step_file;
  string m_track_time_delay_file;

  VEC3D_T m_egy_cmb;

  VEC3D_T m_R_0;
  VEC3D_T m_DE;

  DIRBRBase *m_ebl;
  Table2D *m_optDepthTable;
  double m_tauCutoff;

  RandomNumbers *m_rng;
  MagneticGrid *m_BFieldGrid;
  PairProduction *m_pspace;
  KleinNishina *m_kspace;

  // HDF5 structures:
  struct SecPhoton {
    double egyPrim, egySec, theta, phi, time, thetap, xi, weight;
  };

  struct LeptonData {
    double redshift, egy, px, py, pz, time, rx, ry, rz;
  };

  // Data storage:
  std::vector<SecPhoton> m_secPhotons;
  std::map<std::string, std::vector<LeptonData>> m_leptonData;
};

} // namespace IGCascade

#endif // IGCASCADE_IGCASCADESIM_H
