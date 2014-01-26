/*!
   \file    DIRBR.cpp

    Diffuse InfraRed Background Radiation
 
    \author Vladimir V. Vassiliev \n
            Department of Physics and Astronomy \n
            UCLA \n
            E-Mail: vvv@astro.ucla.edu
	    
    \author Timothy C. Arlen \n
            Department of Physics and Astronomy \n
            UCLA \n
            E-Mail: arlen@astro.ucla.edu

    \date   November 13, 2004

    \version 1.0

    \revision - 5/25/2012 - TCA, rearranged constants, embedded within
    inherited class structure did not add/subtract anything to any
    major class methods.
 
    \note

 */

#ifndef   DIRBRH                  // Header File Guard
#define   DIRBRH
#include "DIRBR.hpp"
#endif                            // DIRBRH


//-------------------------------------------------------------------
// constructors /////////////////////////////////////////////////////
//-------------------------------------------------------------------
//----------------- DIRBR:: DIRBR() --------------------------------- 
/*! Default constructor 
    Sets initialization parameters and default integration steps for 
    redshift and dlambda/lambda

 */
DIRBR::DIRBR( )
{
  K=0;                // initialization of dimension of DIRBR model
  Nq_min=0;           // initialization of dimension of auxiliary arrays

  dlnlambda=0.0002;   // integration step for dln(wavelength)
  q_min=5.0e-4;       // range of definition for integration [q_min, 1]

  //dz=0.001;          // integration step for redshift
  dz = 1.0e-5;

  m_Hubble = 0.7;
  m_OmegaM = 0.3;
  m_OmegaR = 8.4e-5;
  m_OmegaL = 0.7-m_OmegaR;
  double Ho = m_Hubble*100.0/kmToMpc_d;
  // Scaling constant for DIRBR intensity
  nuFnu=(8./3.)*(Ho/CGS_THOM_CS_d/(4.*M_PI))*(meC2_d*meC2_d/ergToTeV_d);
                  // [erg s^-1 cm^-2 sr^-1]
  nuFnu*=1.e+6;   // [erg s^-1 cm^-2 sr^-1] -> [nW m^-2 sr^-1]

  // Scaling constant for wavelength
  kappa=ergCM_hc_d*(ergToTeV_d/meC2_d/meC2_d);   // [cm]
  kappa*=1.e+04;              // [cm] -> [\mum]

  double gamma_=2.5;          // default spectral index for power law SED
  SetDIRBRPowerLaw( gamma_ );

  DIRBR_err=SUCCESS;  // initialization of error status
}
//-------------------------------------------------------------------
/*! Overloaded constructor

    \param StepdLnWavelength - Integration step for dlambda/lambda

    \param StepRedShift      - Integration step for redshift

 */
DIRBR::DIRBR( double StepdLnWavelength, double StepRedShift )
{
  K=0;                // initialization of dimension of DIRBR model
  Nq_min=0;           // initialization of dimension of auxiliary arrays

  dlnlambda=StepdLnWavelength;  // integration step for dln(wavelength)
  q_min=5.0e-4;       // range of definition for integration [q_min, 1]

  dz=StepRedShift;    // integration step for redshift 

  m_OmegaM = 0.3;
  m_OmegaR = 8.4e-5;
  m_OmegaL = 0.7-m_OmegaR;
  m_Hubble = 0.7;
  double Ho = m_Hubble*100.0/kmToMpc_d;
  // Scaling constant for DIRBR intensity
  nuFnu=(8./3.)*(Ho/CGS_THOM_CS_d/(4.*M_PI))*(meC2_d*meC2_d/ergToTeV_d);  
                  // [erg s^-1 cm^-2 sr^-1]
  nuFnu*=1.e+6;   // [erg s^-1 cm^-2 sr^-1] -> [nW m^-2 sr^-1]

  // Scaling constant for wavelength
  kappa=ergCM_hc_d*(ergToTeV_d/meC2_d/meC2_d);   // [cm]
  kappa*=1.e+04;              // [cm] -> [\mum]

  double gamma_=2.5;          // default spectral index for power law SED
  SetDIRBRPowerLaw( gamma_ );

  DIRBR_err=SUCCESS;  // initialization of error status
}
//-------------------------------------------------------------------
// destructor ///////////////////////////////////////////////////////
//-------------------------------------------------------------------
/*! 

 */
DIRBR::~DIRBR()
{

  // Diagnostic print out. 
  /*
  cout<<endl;
  cout<<"sigma240 "<<sigma240<<endl;
  cout<<"Imax240 "<<Imax240<<endl;
  cout<<"Imin240 "<<Imin240<<endl<<endl;

  cout<<"I1216A "<<I1216A<<endl;
  cout<<"Imax1216A "<<Imax1216A<<endl;
  cout<<"Imin1216A "<<Imin1216A<<endl;
  cout<<"IndexUV "<<IndexUV<<endl<<endl;

  cout<<"K "<<K<<endl;
  cout<<"DIRBR_err "<<DIRBR_err<<endl<<endl;

  cout<<"nuFnu: "<<nuFnu<<endl;
  cout<<"kappa: "<<kappa<<endl;
  cout<<"Nq_min: "<<Nq_min<<endl<<endl;

  cout<<"dlnlambda "<<dlnlambda<<endl;
  cout<<"dz "<<dz<<endl<<endl;

  cout<<"gamma "<<gamma<<endl<<endl;

  if( K > 0 ) {
    for (int j=0; j<K; j++) {
      cout<<lambda[j]<<" um"<<I[j]<<" nW/cm^2/sr"<<endl;
    } 
  }
  cout<<endl;
  */

  DeleteMemoryS();
  DeleteMemoryT();
}
//-------------------------------------------------------------------
// accessor Set functions ///////////////////////////////////////////
//-------------------------------------------------------------------
//--------------------- DIRBR::SetS() -------------------------------
/*
! 
   To initialize DIRBR the first and the last elements of wavelength
   array MUST be 240 \mum and 1216 A (0.1216 \mum) respectively. 
   Array of wavelengths must be also ordered from far IR to extreme UV
   and SED array must be positevly defined (>0). 

    \param Kwb - number of wavelength bands [1]

    \param lambda_ - array of wavelength in 10^-6 m [um]

    \param I_ - Spectral Energy Density [nW m^{-2} sr^{-1}]

    \return DIRBR_ERR - error code

    \sa GetDIRBR

 */
DIRBR_ERR DIRBR::SetDIRBR(int Kwb, double * lambda_, double * I_ )
{

  // Check if DIRBR was initialized already and if so deallocate memory
  if ( K > 0 ) DeleteMemoryS();

  // Check if there are at least two points in DIRBR array
  DIRBR_err=ERR_003;
  if ( Kwb < 2 ) return DIRBR_err; 

  // Check for the first and the last wavelength elements in array 
  DIRBR_err=ERR_004;
  if ( lambda_[0] != 240. ) return DIRBR_err;
  DIRBR_err=ERR_005;
  if ( lambda_[Kwb-1] != 0.1216 ) return DIRBR_err;

  // Check for ordering of wavelength array
  DIRBR_err=ERR_006;
  int j;
  for ( j=0; j<Kwb-1; j++ ) {
    if( lambda_[j] < lambda_[j+1] ) return DIRBR_err;
  }

  // Check for positive definition of SED
  DIRBR_err=ERR_007;
  for ( j=0; j<Kwb; j++ ) {
    if( I_[j] <= 0. ) return DIRBR_err;
  }

  // Check for compatible boundary conditions at 240 and 0.1216 \mum
  DIRBR_err=SetFIRAS240(I_[0]);
  if(DIRBR_err != SUCCESS) return DIRBR_err;
  DIRBR_err=Set1216A(I_[Kwb-1]);
  if(DIRBR_err != SUCCESS) return DIRBR_err;

  // Allocate and initialize spline arrays
  K=Kwb;
  lambda = new double [K];
  I = new double [K];
  Lnlambda = new double [K];
  LnI = new double [K];
  h = new double [K];
  m = new double [K];
 
  for (j=0; j<K; j++) {
    lambda[j]=lambda_[j];
    I[j]=I_[j];
    Lnlambda[j]=log(lambda[j]);
    LnI[j]=log(I[j]);
  }

  for (j=0; j<K-1; j++) {
    h[j]=Lnlambda[j]-Lnlambda[j+1];
  }

  // Cubic spline
  {
    double * A;
    double * B;
    A = new double [K];
    B = new double [K];

    A[0]= -1./2.;
    for (j=1; j<K-1; j++) {
      A[j]=-h[j]/(2*h[j]+2*h[j-1]+h[j-1]*A[j-1]);
    }

    // first derivative is continuous with FIRAS fit
    double yp;
    yp=log(GetFIRAS(241.)/GetFIRAS(240.))/log(241./240.);

    B[0]=3*(yp - (LnI[0]-LnI[1])/h[0])/h[0];
    for (j=1; j<K-1; j++) {
      B[j]=A[j]*(
           B[j-1]*h[j-1]/h[j]
          -6*((LnI[j-1]-LnI[j])/h[j-1]-(LnI[j]-LnI[j+1])/h[j])/h[j]);
    }

    // second derivative in UV is zero (power law)
    m[K-1]=0.; 
    for (j=K-1; j>0; j--) m[j-1]=A[j-1]*m[j]+B[j-1];

    // continuous first derivative into Lyman region
    // IndexUV=(LnI[K-1]-LnI[K-2])/h[K-2]+m[K-2]*h[K-2]/6.;

    delete [] A;
    delete [] B;
  }

  DIRBR_err=SUCCESS;
  return DIRBR_err;
}
//-------------------------------------------------------------------
/*! Sets spectral index for power law DIRBR differential photon 
    number density 

    \param gamma  - differential spectral index of power law dn/de

    \sa GetDIRBRPowerLaw(double) 

 */
void DIRBR::SetDIRBRPowerLaw( double gamma_ )
{
  gamma=gamma_;

  // Allocate optical depth arrays if this is the first call 
  if( Nq_min == 0 ) AllocateMemoryT ();

  PLfactor=SEDFactor();

}
//-------------------------------------------------------------------
// accessor Get functions ///////////////////////////////////////////
//-------------------------------------------------------------------
/*! Calculates SED of DIRBR at wavelength l

    \param lambda - wavelength [\mum]

    \return Spectral Energy Density I [nW m^{-2} sr^{-1}]

    \sa SetDIRBR

 */
double DIRBR::GetDIRBR( double l, double z )
/*
  note: z is not used. It is here, only becaue it needs to be the same function
        call as the DIRBR_zEvol class, to be compatible with the 
	airProduction::PropagateEBLRedshift() function.
*/
{
  double LnSED;

  if( K == 0 ) return 0.;

  // extreme UV end lambda < 912 A should be absorbed 
  // by hydrogen ionization in local Universe
  if( l < 0.0912 ) return 0.;

  // Lyman transitions region
  if( l < 0.1216 ) {
    LnSED=LnI[K-1]-IndexUV*log(l/lambda[K-1]);
    return exp(LnSED);
  }

  // far IR
  if( l > 240. ) {
    return GetFIRAS(l);
  }

  //cubic spline in lambda[K-1] < l < lambda[0] region
  int j=1;
  while( (lambda[j-1]-l)*(l-lambda[j]) < 0. ) {
    j++;
  }

  double x=log(l/lambda[j])/h[j-1];
  LnSED=(LnI[j  ]-m[j  ]*h[j-1]*h[j-1]/6.)*(1-x)+
        (LnI[j-1]-m[j-1]*h[j-1]*h[j-1]/6.)*   x +
        m[j  ]*h[j-1]*h[j-1]/6. *(1-x)*(1-x)*(1-x)+
        m[j-1]*h[j-1]*h[j-1]/6. *   x *   x *   x;

  return exp(LnSED);
}
//-------------------------------------------------------------------
/*! Calculates power law SED of DIRBR with spectral 
    index 2-gamma.

    \param lambda - wavelength [\mum]

    \return Spectral Energy Density I [nW m^{-2} sr^{-1}]

    \sa SetDIRBRPowerLaw(double) 

 */
double DIRBR::GetDIRBRPowerLaw( double lambda )
{
  return nuFnu*pow(lambda/kappa, gamma-2.);
}
//-------------------------------------------------------------------
/*! Calculates energy in the wavelength interval

    \param lambda1 - wavelength bound 1 [\mum]
    \param lambda2 - wavelength bound 2 [\mum]

    \return Energy [nW m^{-2} sr^{-1}]

    \sa SetDIRBR

 */
double DIRBR::GetEnergy( double lambda1, double lambda2 )
{
  // sanity checks 
  double lambdaLower=lambda1;
  double lambdaUpper=lambda2;
  if(lambdaUpper < lambdaLower ) {
    lambdaUpper=lambda1;
    lambdaLower=lambda2;
  }
  if(lambdaUpper < 0.0912) return 0.;
  if(lambdaLower < 0.0912) lambdaLower = 0.0910;


  double SED=0.0; double Energy;
  double q=exp(dlnlambda);

  double l=lambdaLower;
  Energy=GetDIRBR(l)/2.;

  for (l=lambdaLower; l<lambdaUpper ; l*=q) {
    SED=GetDIRBR(l);
    Energy+=SED;
  }
  Energy-=SED/2.;
  Energy*=dlnlambda;

  return Energy;
}
//-------------------------------------------------------------------
/*! CMBR

    \param lambda - wavelength [\mum]

    \return Spectral Energy Density I [nW m^{-2} sr^{-1}]

    \sa 

 */
double DIRBR::GetCMBR( double& lambda )
{

  double Tcmb=2.728*eV_K_B_d;         // T = 2.728 +/- 0.004 K
  double e=1.2415/lambda;       // \mum -> eV

  return 5.042e+16*e*e*e*e/(exp(e/Tcmb)-1.);
}
//-------------------------------------------------------------------
// Methods //////////////////////////////////////////////////////////
//-------------------------------------------------------------------
/*! This routine tests all direct constraints on DIRBR SED 
    It excludes detections at 0.3, 0.555, and 0.814 reported by 
    Bernstein et al. 2002, ApJ, 571, 56 because allowed error
    bars overlap all possible values limited by the other measurements. 

    \return DIRBR_ERR - error code reflecting violating constraints

 */
DIRBR_ERR DIRBR::TestDIRBRlimits( void )
{
  // 95%CL DIRBE residuals 
  // M.G. Hauser, et al. 1998, ApJ, 508, 25
  DIRBR_err=ERR_008;
  if ( GetDIRBR(1.25)  >  75. ) return DIRBR_err; 
  if ( GetDIRBR( 2.2)  >  39. ) return DIRBR_err; 
  if ( GetDIRBR( 3.5)  >  23. ) return DIRBR_err; 
  if ( GetDIRBR( 4.9)  >  41. ) return DIRBR_err; 
  if ( GetDIRBR( 12.)  > 468. ) return DIRBR_err; 
  if ( GetDIRBR( 25.)  > 504. ) return DIRBR_err; 
  if ( GetDIRBR( 60.)  >  75. ) return DIRBR_err; 
  if ( GetDIRBR(100.)  >  34. ) return DIRBR_err; 
  if ( GetDIRBR(140.)  >  43. ) return DIRBR_err; 
  if ( GetDIRBR(240.)  >  20. ) return DIRBR_err; 

  // 240 \mum Hauser et al. 1998, ApJ,508, 25 
  // Detection: 13+-2 (FIRAS calibration)
  // Detection: 14+-3 (DIRBE calibration)
  // 2 sigma limits:
  //
  // NOTE: (TCA) - changed lower limit from 8 to 3 nW m^-2 sr^-1, 
  //   from reference Berta et al. A&A, 2010, 518, L30
  // 
  // FURTHERMORE: I am removing the DIRBR limit errors at 140.0,
  // 100.0, and 70.0, for basically the same reason (no longer
  // considered robust).
  DIRBR_err=ERR_009;
  if ( GetDIRBR(240.)  >  20. ) return DIRBR_err; 
  if ( GetDIRBR(240.)  <  3.  ) return DIRBR_err;

  // 140 \mum Hauser et al. 1998, ApJ,508, 25 
  // Detection: 15+-6 (FIRAS calibration)
  // Detection: 25+-7 (DIRBE calibration)
  // 2 sigma limits:
  //DIRBR_err=ERR_010;
  //if ( GetDIRBR(140.)  >  39. ) return DIRBR_err; 
  //if ( GetDIRBR(140.)  <  3.  ) return DIRBR_err; 

  // 100 \mum Lagache et al. 2000, A&A, 354, 247 
  // Detection: 23+-6
  // 2 sigma limits:
  //DIRBR_err=ERR_011;
  //if ( GetDIRBR(100.)  >  35. ) return DIRBR_err; 
  //if ( GetDIRBR(100.)  <  11. ) return DIRBR_err; 

  // 70 \mum Frayer, et al. 2006, ApJ, 647, L9
  // Lower limits from galaxy counts using GOODS-North instrument.
  // 2 sigma limits:
  //DIRBR_err=ERR_012;
  //if( GetDIRBR(70.) < 3.6 ) return DIRBR_err;

  // P.B. Hacking \& B.T. Soifer, 1991, ApJ, 367, L49
  // Lower limits from galaxy counts derived from
  // IRAS observations. 100 \mum (>3.89) and 25 \mum (>1.02)
  // 2 sigma limits are now obsolete because of the Lagache 
  // detection, and SPITZER resolved galaxy counts.
  DIRBR_err=ERR_013;
  if ( GetDIRBR(60.)  <  1.89 ) return DIRBR_err;  
 
  // Papovich et al. 2004, ApJS, 154, 70
  // Lower limits from galaxy counts derived from
  // SPITZER surveys at 24 \mum. Resolved 1.9+-0.6, 
  // but extrapolation to fainter flux suggest 2.7+1.1/-0.7   
  DIRBR_err=ERR_014; 
  if ( GetDIRBR(24.)  <  1.3 ) return DIRBR_err; 

  // Elbaz et al. 2002, A&A, 384, 848  :  2.4+/-0.5 at 15 \mum
  // Metcalfe et al. 2003, A&A, 407, 791  : 2.7+/-0.62 at 15 \mum
  // ISO (Infrared Space Observatory) galaxy counts
  DIRBR_err=ERR_015; 
  if ( GetDIRBR(15.)  <  1.4 ) return DIRBR_err; 

  // Madau & Pozzetti 2000, MNRAS, 312, L9  
  // HST galaxy number counts
  // 2.2       7.9+2.0/-1.2
  // 1.6       9.0+2.6/-1.7
  // 1.1       9.7+3.0/-1.9
  // 0.81      8.0+1.6/-0.9
  // 0.67      6.7+1.3/-0.9
  // 0.45      4.6+0.7/-0.5
  // 0.36      2.9+0.6/-0.4
  DIRBR_err=ERR_016; 
  if ( GetDIRBR(2.2)  <  5.5 ) return DIRBR_err; 
  if ( GetDIRBR(1.6)  <  5.6 ) return DIRBR_err; 
  if ( GetDIRBR(1.1)  <  5.9 ) return DIRBR_err; 
  if ( GetDIRBR(0.81) <  6.2 ) return DIRBR_err; 
  if ( GetDIRBR(0.67) <  4.9 ) return DIRBR_err; 
  if ( GetDIRBR(0.45) <  3.6 ) return DIRBR_err; 
  //if ( GetDIRBR(0.36) <  2.1 ) return DIRBR_err;  //Changed-TCA June 10, 2009
  if ( GetDIRBR(0.36) <  1.8 ) return DIRBR_err; 

  // C. Armand et al. 1994, A&A, 284, 12
  // 2000 A  galaxy counts estimate 95%CL lower limit of 0.79
  DIRBR_err=ERR_017; 
  if ( GetDIRBR(0.2)  <  0.79 ) return DIRBR_err; 

  // C. Kevin Xu et al. 2004, astro-ph/0411317
  // GALEX galaxy counts at 1530 A and 2310 A
  // 0.153    1.03+/-0.15 
  // 0.231    2.25+/-0.32 
  DIRBR_err=ERR_018; 
  if ( GetDIRBR(0.231)  <  1.61 ) return DIRBR_err; 
  if ( GetDIRBR(0.153)  <  0.73 ) return DIRBR_err; 
 
  // Gardner et al. 2000 ApJ, 542, L79
  // UV galaxy counts from STIS observations of 
  // the Hubble Deep Fields (caution ! low statistics)
  // NUV; 2365 A    3.6+0.7/-0.5
  // FUV; 1595 A    2.9+0.6/-0.4 (lower limit)
  DIRBR_err=ERR_019; 
  // Disable due to low statistic and inconsistency with 
  // GALEX result and earlier Armand et al. result 
  /*
  if ( GetDIRBR(0.2365)  <  2.6  ) return DIRBR_err; 
  if ( GetDIRBR(0.1595)  <  2.1  ) return DIRBR_err; 
  */

  // Bowyer S., 1991, ARA&A, 29, 59
  // All-sky photometry upper limits as reported
  // in  L. Pozzetti et al. 1998, MNRAS, 298, 1133
  // 0.165    5.56 
  // see also astro-ph/0004147 for 0.145-0.19 \mum 
  // conservative upper limit from STIS on HST which 
  // is 9.94+/-2.04 
  DIRBR_err=ERR_020; 
  if ( GetDIRBR(0.165)  >  5.56 ) return DIRBR_err; 

  // M. Maucherat-Joubert et al. 1980, A&A, 88, 323
  // Upper limit (95%CL) left as a residual after all
  // substractions. The data from UV survey performed
  // by the ELZ photometer on board the D2B satellite
  // are used for this analysis.
  // 0.220   6.9
  DIRBR_err=ERR_021; 
  if ( GetDIRBR(0.220)  >  6.90 ) return DIRBR_err; 

  // G.N. Toller 1983, ApJ, 266, L79
  // Upper limit is established using the
  // photopolarimeter on board the Pioneer 10
  // spacecraft.
  // 0.440  19.8
  DIRBR_err=ERR_022; 
  if ( GetDIRBR(0.440)  >  19.8 ) return DIRBR_err; 
 
  // R. R. Dube et al. 1979, ApJ, 232, 333
  // Upper limit obtained from photoelectric
  // measurements of the separate sources
  // contributing to the night sky brightness
  // at 0.5115 $\mu$m.
  // 0.5115  26.0
  DIRBR_err=ERR_023; 
  if ( GetDIRBR(0.5115)  >  26.0 ) return DIRBR_err; 

  DIRBR_err=SUCCESS;
  return DIRBR_err;
}
//-------------------------------------------------------------------


double DIRBR::OpticalDepth( double& E, double& RedShift )
/*! Calculates optical depth for a VHE photon traveling from 
    a source at redshift z. No redshift evolution of DIRBR is 
    assumed except energy and density rescaling with appropriate
    (1+z) factors.

    \param E - detected by observer photon energy [TeV]

    \param z - redshift of a source [1]

    \return  optical depth  [1] 

 */
{
  if ( K == 0 ) return 0.; // DIRBR has not been initialized

  // Allocate optical depth arrays if this is the first call 
  if( Nq_min == 0 ) AllocateMemoryT ();

  // find optical depth

  // Redshift integral
  double z=0.;
  double Tau=0.;

  //std::ofstream tau_integral("tau_integral_DIRBR.txt");

  for (z=0.; z < RedShift ; z+=dz) { 

    double z1=(1.+z)*1.;
    double z2=(1.+z)*z1;
    double z3=(1.+z)*z2;
    double z4=(1.+z)*z3;

    double l_max=kappa*E*z2;
    double W=sqrt(m_OmegaM*z3+m_OmegaL+(1.-m_OmegaM-m_OmegaL)*z2);

    // Wavelength integral
    double tau=0.;
    double l,dbr=1.;
    for (int i=Nq_min-1; i>=0 && dbr != 0.; i--) {
      l=l_max*q[i];

      // Find energy density 
      dbr=GetDIRBR(l);
      if(l > 100.) dbr+=GetCMBR(l);

      // Test optical depth integration with power law SED
      // dbr=GetDIRBRPowerLaw(l);

      tau+=dbr*F3[i];
    }
    tau*=dlnlambda/nuFnu;
    tau*=E;
    tau*=z4/W;

    if (dz < RedShift-z) {
      Tau+=tau*dz;
    } else {
       Tau+=tau*(RedShift-z);
    }

    //tau_integral<<"tau: "<<Tau<<" z: "<<z<<std::endl;
    
  }

  //tau_integral.close();

  return Tau;
}
//-------------------------------------------------------------------
/*! Calculates optical depth for a VHE photon traveling from 
    a source at redshift z for power law DIRBR SED. No redshift 
    evolution of DIRBR is assumed except energy and density 
    rescaling with appropriate (1+z) factors.

    \param E - detected by observer photon energy [TeV]

    \param z - redshift of a source [1]

    \return  optical depth  [1] 

 */
double DIRBR::OpticalDepthPowerLaw( double E, double RedShift )
{

  return pow(E, gamma-1.)*RedShiftFactor(RedShift)*PLfactor;

}
//-------------------------------------------------------------------
/*! Calculates distance to a source at redshift z. No redshift
    evolution of DIRBR is assumed.

    \param z - redshift of a source [1]

    \return  distance  [1] 

 */
double DIRBR::Distance( double& RedShift )
{
  // Redshift integral
  double z=0.;
  double d=0.;

  for (z=0.; z < RedShift ; z+=dz) { 

    double z1=(1.+z)*1.;
    double z2=(1.+z)*z1;
    double z3=(1.+z)*z2;
    double z4=(1.+z)*z3;

    double W=sqrt(m_OmegaR*z4+m_OmegaM*z3+m_OmegaL+(1.-m_OmegaR-m_OmegaM-m_OmegaL)*z2);

    if (dz < RedShift-z) {
      d+=1.0/W*dz;
    } else {
      d+=1.0/W*(RedShift-z);
    }
  }

  double Ho = m_Hubble*100.0/kmToMpc_d;
  d*=CGS_C_d/Ho;
  return d;
}
//-------------------------------------------------------------------
// Private Methods //////////////////////////////////////////////////
//-------------------------------------------------------------------
/*! Allocat memory and initialize internal arrays 
    for calculation of optical depth 
 */
void DIRBR::AllocateMemoryT()
{
    Nq_min=(int)(-log(q_min)/dlnlambda);
    q = new double [Nq_min];
    f = new double [Nq_min];
    F3= new double [Nq_min];

    double q0=exp(dlnlambda);
    double y=exp(-(Nq_min-1)*dlnlambda);
    
    for (int i=0; i<Nq_min-1; i++) {
      q[i]=y;
      f[i]=y*(
	      (1.+y-y*y/2.)*log((1.+sqrt(1.-y))/(1.-sqrt(1.-y)))
	      -(1.+y)*sqrt(1.-y)
	      );
      y*=q0;
    }
    q[Nq_min-1]=1.;
    f[Nq_min-1]=0.;

    F3[Nq_min-1]=0.;
    for (int i=Nq_min-2; i>=0; i--) {
      F3[i]=F3[i+1]+dlnlambda*f[i]/q[i]/q[i];
    }
    for (int i=0; i<Nq_min; i++) {
      F3[i]*=2*q[i]*q[i]*q[i];
    }
}
//-------------------------------------------------------------------
/*! Delete allocated memory for auxiliary optical depth arrays

 */
void DIRBR::DeleteMemoryT()
{

  if( Nq_min > 0 ) {
    delete [] q;
    delete [] f;
    delete [] F3;
  }
  Nq_min=0;

}
//-------------------------------------------------------------------
/*! Delete allocated memory for internal spline arrays

 */
void DIRBR::DeleteMemoryS()
{
  if( K > 0 ) {
    delete [] lambda;
    delete [] I;
    delete [] Lnlambda;
    delete [] LnI;
    delete [] h;
    delete [] m;
  }
  K=0;
}
//-------------------------------------------------------------------
/*! Finds energy and impact anlgle integral for power law DIRBR SED with 
    differential photon density spectral index - gamma 

    \return factor  [1]

    \sa GetDIRBRPowerLaw(double), SetDIRBRPowerLaw(double)

 */
double DIRBR::SEDFactor( void )
{
  if ( K == 0 ) return 0.; // DIRBR has not been initialized

  // Allocate arrays if this is the first call for a given DIRBR 
  if( Nq_min == 0 ) AllocateMemoryT ();

  double y=0.;
  for (int i=0; i<Nq_min; i++) {
    y+=f[i]*pow(q[i],gamma)/q[i];      
  }
  y*=dlnlambda;

  y*=2./(gamma+1.);   // contribution from integration over angles

  return y;

}
//-------------------------------------------------------------------
/*! Finds redshift integral for power law DIRBR SED with 
    differential photon density spectral index - gamma

    \param RedShift - red shift [1]

    \return factor  [1]

    \sa GetDIRBRPowerLaw(double), SetDIRBRPowerLaw(double)
          
 */
double DIRBR::RedShiftFactor( double RedShift )
{

  // Redshift integral
  double z;
  double factor;
  double z1,z2,z3,zg;
  double W;

  factor=0.;
  for (z=0.; z < RedShift ; z+=dz) { 

    z1=(1.+z)*1.;
    z2=(1.+z)*z1;
    z3=(1.+z)*z2;
    zg=pow(1.+z, 2*gamma);

    W=sqrt(m_OmegaM*z3+m_OmegaL+(1.-m_OmegaM-m_OmegaL)*z2);

    if (dz < RedShift-z) {
      factor+=zg*dz/W;
    } else {
       factor+=zg*(RedShift-z)/W;
    }
  }

  return factor;
}
//-------------------------------------------------------------------
// accessor Set functions ///////////////////////////////////////////
//-------------------------------------------------------------------
//-------------------------------------------------------------------
/*! Sets SED at 1216 A. 

    The bounds in this routine are based on power law 
    extrapolation (I~(lambda)^1.33) from upper and low bounds 
    in near and middle UV which produce 0.38 and 4.0 estimates.

 */
DIRBR_ERR DIRBR::Set1216A( double I )
{
  Imax1216A=6.0;
  Imin1216A=0.2;
  IndexUV=-2.5;

  if(I < Imin1216A || I > Imax1216A ) { 
    I1216A=1.;
    DIRBR_err=ERR_002;
    return DIRBR_err;
  }

  I1216A=I;
  DIRBR_err=SUCCESS;
  return DIRBR_err;
}
//-------------------------------------------------------------------
/*! Sets DIRBR SED amplitude at 240 \mum 

    \param I - DIRBR SED at 240 \mum [nW m^{-2} sr^{-1}]

    \return 'false' if incompatible with FIRAS 2sigma bounds
            3.01574 < I <  21.8873
            D.J. Fixsen et al. (1998), ApJ, 508, 123
            otherwise returns 'true'

 */
DIRBR_ERR DIRBR::SetFIRAS240( double I )
{
  int k;
  double Io;
  double lb=2.;
  double rb=-2.;
  double lambda=240.;
  sigma240=lb;
  Imin240=GetFIRAS(lambda);
  sigma240=rb;
  Imax240=GetFIRAS(lambda);
  
  if(I < Imin240 || I > Imax240 ) { 
    sigma240=0.;
    DIRBR_err=ERR_001;
    return DIRBR_err;
  }
  
  for (k = 0; k < 32; k++) {
    sigma240=(lb+rb)/2.;
    Io=GetFIRAS(lambda);
    if(Io < I ) lb=sigma240;
    if(Io > I ) rb=sigma240;
  }
  
  DIRBR_err=SUCCESS;
  return DIRBR_err;
}
//-------------------------------------------------------------------
// accessor Get functions ///////////////////////////////////////////
//-------------------------------------------------------------------
//-------------------------------------------------------------------
/*! COBE FIRAS result in 125 - 2000 \mum region 
    D.J. Fixsen et al. (1998), ApJ, 508, 123

    \param lambda - wavelength [\mum]

    \return Spectral Energy Density I [nW m^{-2} sr^{-1}]

    \sa SetFIRAS240( double & )

 */
double DIRBR::GetFIRAS( double lambda )
{
  // applicability condition
  if(lambda < 125.) return 0.; 
  if(lambda > 2000.) return 0.; 

  double T=(18.5-sigma240*1.2)*eV_K_B_d;
  double Index=4.64+sigma240*0.12;  
  double Amplitude=1.+sigma240*0.30769;

  double e=1.2415/lambda;       // \mum -> eV

  double p=1./(exp(e/T)-1.);

  return Amplitude*1.0886e+13*pow(e,Index)*p;
}
