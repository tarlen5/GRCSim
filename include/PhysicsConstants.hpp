/*!
---------------------------------------------------------------------------
  \file     PhysicsConstants.hpp

	    Brief constants class which go under the namespace PhysicsConstants

  \author   Timothy C. Arlen                          \n
            Department of Physics and Astronomy       \n
            UCLA                                      \n
	    arlen@astro.ucla.edu                      \n

  \date     March 16, 2010                            \n
---------------------------------------------------------------------------
*/

#include <qd/dd_real.h>
#include <Vec3D.hpp>

#ifndef PHYSICSCONSTANTS_H
#define PHYSICSCONSTANTS_H

namespace PhysConst {

  // Physical Constants:
  const VEC3D_T CGS_C     = "2.99792458E+10";      // [cm/s]
  const VEC3D_T MPC_TO_CM = "3.086E+24";           // [cm Mpc^-1]  
  const VEC3D_T SI_MELEC  = "9.109382E-28";      // [g]
  const VEC3D_T nJ_TO_eV  = "6.241506E+09";              
  const VEC3D_T J_TO_eV   = nJ_TO_eV*1.0E+9;             
  const VEC3D_T eV_MELEC  = SI_MELEC*CGS_C*CGS_C*1.E-7*J_TO_eV; // [eV]
  const VEC3D_T TeV       = "1.0E+12";
  const VEC3D_T CM_TO_PC  = "3.086E+18";        // [cm pc^-1]
  const VEC3D_T eV_K_B    = "8.61734315E-5";   // [eV/K] Boltzmann const
  // [esu]
  const VEC3D_T ESU_ELEC_CHARGE = "4.80320440498320182E-10"; 
  // [eV*s] Planck cons
  const VEC3D_T eVSEC_PLANCK = "4.1356674335E-15";
  // planck const * speed of light = hc [eV cm]
  const VEC3D_T eVCM_HC = eVSEC_PLANCK*CGS_C; 
  // [cm] Classical elec rad:
  const VEC3D_T CGS_ELEC_RAD = 
    ESU_ELEC_CHARGE*ESU_ELEC_CHARGE/SI_MELEC/CGS_C/CGS_C;   
  // Thomson c-s [cm^2]
  const VEC3D_T CGS_THOM_CS = 8.0/3.0*VEC3D_PI*CGS_ELEC_RAD*CGS_ELEC_RAD;


  // LCDM FRW Cosmology
  const VEC3D_T HUB_CONST = "2.2683e-18";
  const VEC3D_T CGS_HUBRAD = CGS_C/HUB_CONST;   // [cm] Hubble Radius
  const VEC3D_T OMEGA_R = "8.4e-5";             // Omega radiation
  const VEC3D_T OMEGA_M = "0.3";                // Omega matter   
  const VEC3D_T OMEGA_L = 0.7 - OMEGA_R;        // Omega lambda
  const VEC3D_T OMEGA_0 = OMEGA_M+OMEGA_R+OMEGA_L;
  
  
}

#endif // PHYSICSCONSTANTS_H
