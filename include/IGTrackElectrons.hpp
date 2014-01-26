/*!
-------------------------------------------------------------------------------
    \file   IGTrackElectrons.hpp
    
    Header file for IGTrackElectrons.cpp
  
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
  
  class IGTrackElectrons
  {
  public:
    
    IGTrackElectrons(const string& s_egy, const string& s_Bmag, const string& s_ze,
		 const string& s_cellsize, const string& s_iterations, 
		 const string& s_eblmodel_file, const string& s_file_count);

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
    void DefineStaticVariables(const string& static_var_file);
    void DefineRp(void);
    void InitializeEBL(const string& ebl_model_file);

    void CreateDirectGamma(RelParticle& GammaPhoton);
    bool PropagateDirectPhoton(RelParticle& GammaPhoton, RelParticle* Electron, 
			 RelParticle* Positron);
    void DefineLeptons(RelParticle& GammaPhoton, RelParticle& EBLPhoton, 
		       RelParticle* Electron, RelParticle* Positron);
    
    bool PropagateSecondaryPhoton(RelParticle& GammaPhoton, 
				  RelParticle* Electron, 
				  RelParticle* Positron);
    double GetMinZ(const Vec4D& gam_ph_p4, const Vec4D& gam_ph_r4, 
		   const VEC3D_T& gam_ph_z);
    void PropagateLepton(RelParticle* Lepton,RelParticle& GammaPhoton);
    void SaveToLowEnergyFile(RelParticle& Particle);
    void SaveDirectPhoton(RelParticle& Particle, std::ofstream& photon_list);
    void SaveLepton(RelParticle* Lepton);

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
    
    // From StaticVariables:
    VEC3D_T m_egy_cmb;
    VEC3D_T m_egy_gamma_min;
    VEC3D_T m_egy_lepton_min;
    bool    m_LOCK;

    VEC3D_T m_R_0;

    DIRBRBase*     m_ebl;

    RandomNumbers*  m_rng;
    MagneticGrid*   m_BFieldGrid;
    PairProduction* m_pspace;
    KleinNishina*   m_kspace;
    
  };
  
}

#endif // IGCASCADE_IGCASCADESIM_H
