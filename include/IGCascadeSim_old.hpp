/*!
-------------------------------------------------------------------------------
    \file   IGCascadeSim.hpp
    
    Header file for IGCascadeSim.cpp
  
    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
	       arlen@astro.ucla.edu                  \n

    \date      May 20, 2012                          \n
    
    \note 
-------------------------------------------------------------------------------
*/

#ifndef IGCASCADE_IGCASCADESIM_H
#define IGCASCADE_IGCASCADESIM_H

#include "KleinNishina.hpp"
#include "PairProduction.hpp"
#include "RelParticle.hpp"
#include "Vec4D.hpp"
#include "Vec3D.hpp"
#include "RandomNumbers.hpp"
#include "GalacticGrid.hpp"
#include "PhysicsConstants.hpp"
#include "DIRBR.hpp"
//#include "Table.hpp"

// ROOT INCLUDES:
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stack>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <stdio.h>

//#include <qd/dd_real.h>

using namespace std;

namespace IGCascade
{
  
  class IGCascadeSim
  {
  public:
    
    IGCascadeSim(const string& s_egy, const string& s_Bmag, const string& s_ze,
		 const string& s_cellsize, const string& s_iterations, 
		 const string& s_eblmodel_file, const string& s_file_count);
    inline VEC3D_T GetMinLepEgy() { return m_egy_lepton_min; }

    void CreateDirectGamma(RelParticle& GammaPhoton);
    bool PropagateDirectPhoton(RelParticle& GammaPhoton, RelParticle* Electron, 
			       RelParticle* Positron);
    
    void RunCascade(void);
    
  private:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Private Functions   ///////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    string DefineCascadeFile(const string& s_eblmodel,const string& s_egy,
		      const string& s_Bmag, const string& s_ze, 
		      const string& s_cellsize, const string& s_file_num);
    string DefineMFfile(const string& static_var_file, const string& s_Bmag,
			const string& s_cellsize, const string& s_ze);
    string DefineLowEgyFile(const string& s_eblmodel,const string& s_egy,
		      const string& s_Bmag, const string& s_ze, 
		      const string& s_cellsize, const string& s_file_num);
    string DefineTrackLeptonFile(const string& s_eblmodel, const string& s_egy, 
			const string& s_Bmag, const string& s_ze, 
			const string& s_cellsize, const string& s_file_num);
    string DefineTrackTimeDelayFile(const string& s_eblmodel,
				    const string& s_egy,
				    const string& s_Bmag, const string& s_ze, 
			  const string& s_cellsize, const string& s_file_num);
    void DefineStaticVariables(const string& static_var_file);
    void DefineRp(void);
    void InitializeEBL(const string& ebl_model_file);
    void DefineOptDepthTable(const string& s_eblmodel, const string& s_ze);
    //void DefineOpticalDepthTable(void);
    void DefineLeptons(RelParticle& GammaPhoton, RelParticle& EBLPhoton, 
		       RelParticle* Electron, RelParticle* Positron);
    void CreateLeptonTree(RelParticle* Lepton);
    void WriteLeptonToFile(RelParticle* Lepton);
    double GetOpticalDepthVal(double egy, double z);
    void PropagatePhotonToObserver(RelParticle& GammaPhoton);
    bool PropagateSecondaryPhoton(RelParticle& GammaPhoton, 
				  RelParticle* Electron, 
				  RelParticle* Positron);
    double GetMinZ(const Vec4D& gam_ph_p4, const Vec4D& gam_ph_r4, 
		   const VEC3D_T& gam_ph_z);
    void PropagateLepton(RelParticle* Lepton,RelParticle& GammaPhoton);
    void StoreSecPhoton(Vec4D& gam_ph_p4, Vec4D& gam_ph_r4, 
			const double& weight);
    void SaveToLowEnergyFile(RelParticle& Particle);
    void SaveDirectPhoton(RelParticle& Particle, std::ofstream& photon_list);
    void SaveLepton(RelParticle* Lepton);
    void SaveToTrackTimeDelayFile(RelParticle& Particle);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Private Data Members///////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEC3D_T m_egy_cascade;    // TeV
    VEC3D_T m_bmag;
    VEC3D_T m_ze;    
    string  m_cellsize;
    int     m_iterations;
    
    string  m_cascade_file;
    string  m_low_egy_file;
    string  m_save_lepton_file;
    string  m_track_lepton_step_file;
    string  m_track_time_delay_file;
    
    // From StaticVariables:
    VEC3D_T m_egy_cmb;
    VEC3D_T m_egy_gamma_min;
    VEC3D_T m_egy_lepton_min;
    bool    m_LOCK;

    VEC3D_T m_R_0;
    VEC3D_T m_DE;

    DIRBRBase*     m_ebl;
    Table2D*       m_optDepthTable;
    double         m_tauCutoff;

    RandomNumbers*  m_rng;
    MagneticGrid*   m_BFieldGrid;
    PairProduction* m_pspace;
    KleinNishina*   m_kspace;

    // TTree for secondary photons and for tracking leptons:
    TTree*   m_secPhotonTree;
    Double_t m_egyPrim, m_egySec, m_theta, m_phi, m_time, m_thetap, m_xi, 
      m_weight;
    
  };
  
}

#endif // IGCASCADE_IGCASCADESIM_H
