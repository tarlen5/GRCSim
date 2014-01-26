/*! 
  \file DIRBRBase.hpp
        DIRBR base class header file
  
  \author   Tim Arlen                  \n
            UCLA                       \n
	    arlen@physics.ucla.edu     \n

  \date     July 27, 2010
  \version  0.0
  \note
*/

#ifndef DIRBRBASE_H
#define DIRBRBASE_H

#include <cmath>

enum DIRBR_ERR {SUCCESS,
                ERR_001,   //!< outside of FIRAS bounds at 240 \mum
                ERR_002,   //!< outside of bounds at 1216 A
                ERR_003,   //!< < 2 wavelengths in DIRBR initialization
                ERR_004,   //!< first wavelength is NOT 240 \mum
                ERR_005,   //!< last wavelength is NOT 1216A 
                ERR_006,   //!< wavelength array is not ordered 
                ERR_007,   //!< negative elements in SED initialization
                ERR_008,   //!< above DIRBE residuals
                ERR_009,   //!< outside of bounds at 240\mum DIRBE detection
                ERR_010,   //!< outside of bounds at 140\mum DIRBE detection
                ERR_011,   //!< outside of bounds at 100\mum Lagache detection
                ERR_012,   //!< below GOODS-North Lower Limit at 70 \mum
		ERR_013,   //!< below IRAS galaxy counts at 60 \mum
                ERR_014,   //!< below SPITZER galaxy counts at 24 \mum
                ERR_015,   //!< below ISO galaxy counts at 15 \mum
                ERR_016,   //!< below HST galaxy counts at 0.36 - 2.2 \mum
                ERR_017,   //!< below galaxy counts at 2000 A
                ERR_018,   //!< below GALEX galaxy counts at 1530 A and 2310 A
                ERR_019,   //!< below STIS galaxy counts at 2365A & 1595A
                ERR_020,   //!< above all-sky photometry UL at 1650 A
                ERR_021,   //!< above D2B satellite photometry at 2200 A 
                ERR_022,   //!< above Pioneer 10 result at 4400 A 
                ERR_023,   //!< above UL from night sky brightness at 5115 A
                ERR_024,    //!< 
              };


/////////////////////////////////////////////////////////
/// Physics Constants
/////////////////////////////////////////////////////////
// Physical Constants:
const double CGS_C_d     = 2.99792458E+10;      // [cm/s]
//const double MPC_TO_CM = 3.086E+24;           // [cm Mpc^-1]  
const double SI_MELEC_d  = 9.109382E-28;      // [g]
const double J_TO_eV_d   = 6.241506E+09*1.0E+9;
//const double nJ_TO_eV  = 6.241506E+09;   
//const double J_TO_eV   = nJ_TO_eV*1.0E+9;             
const double eV_MELEC_d  = SI_MELEC_d*CGS_C_d*CGS_C_d*1.E-7*J_TO_eV_d; // [eV]
const double eV_K_B_d    = 8.61734315E-5;   // [eV/K] Boltzmann const
const double ESU_ELEC_CHARGE_d = 4.80320440498320182E-10; 
const double eVSEC_PLANCK_d = 4.1356674335E-15;
// planck const * speed of light = hc [eV cm]
const double eVCM_HC_d = eVSEC_PLANCK_d*CGS_C_d; 
const double CGS_ELEC_RAD_d = 
  ESU_ELEC_CHARGE_d*ESU_ELEC_CHARGE_d/SI_MELEC_d/CGS_C_d/CGS_C_d;   
// Thomson c-s [cm^2]
const double CGS_THOM_CS_d = 8.0/3.0*M_PI*CGS_ELEC_RAD_d*CGS_ELEC_RAD_d;
const double kmToMpc_d = 3.08568e19;
const double ergToTeV_d=1.6;    // erg; 1 TeV = 1.6 erg
const double meC2_d=0.8175982e-06;  // erg electron mass
static const double ergCM_hc_d=1.9864085e-16;    // erg cm

//------------------------------------------
// DIRBRBase - Abstract Base Class
//------------------------------------------
class DIRBRBase
{
public:
  
  DIRBRBase()
  {
    // nothing to see here.
  }

  // -------------------------------------
  // Public Member Functions:
  // -------------------------------------
  // Pure virtual functions:
  virtual double GetDIRBR(double, double) = 0;
  virtual double GetCMBR(double&) = 0;
  virtual double OpticalDepth(double&, double&) = 0;
  virtual double Distance (double&) = 0;
  virtual DIRBR_ERR SetDIRBR(int,double*,double*) = 0;     //!< set DIRBR model 
  virtual DIRBR_ERR TestDIRBRlimits(void) = 0; //!< test existing direct limits
  virtual inline void SetCosmology(const double OmegaM, const double OmegaR, 
				   const double OmegaL, const double Hubble_h) {
    m_OmegaM = OmegaM;
    m_OmegaR = OmegaR;
    m_OmegaL = OmegaL;
    m_Hubble = Hubble_h;
  }

  virtual ~DIRBRBase() { /*virtual destructor*/ }

protected:
  
  double m_OmegaM;
  double m_OmegaR;
  double m_OmegaL;
  double m_Hubble;

};


#endif // DIRBRBASE_H
