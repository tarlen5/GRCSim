/*!
-------------------------------------------------------------------------------
    \file   DIRBR_Redshift.cpp

     DIRBR Redshift Evolution class implementation file.
  
    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
	       arlen@astro.ucla.edu                  \n

    \date      July 24, 2010
  
    \version:  1.0

    \revision: 

    \note:

-------------------------------------------------------------------------------
*/

#include "DIRBR_Redshift.hpp"

DIRBR_Redshift::DIRBR_Redshift(std::string filename)
/*
  IMPORTANT: Format of filename:
  Filename must be a 2D table of format:
  lambda [mum] Intensity(z0) [nW m^-2 sr^-1] Intensity(z1) Intensity(z2) ...

  where the first line is:
  -1 z0 z1 z2...

*/
{
  
  SEDTable = new Table2D(filename);
  
  int num_rows = SEDTable->GetNRows();
  int num_cols = SEDTable->GetNCols();
  m_lambda_min = SEDTable->m_table[1][0];
  m_lambda_max = SEDTable->m_table[num_rows-1][0];
  m_zmin = SEDTable->m_table[0][1];
  m_zmax = SEDTable->m_table[0][num_cols-1];
  
  InitializeIntegral();
  m_dz = 1.0e-5;

  m_OmegaM = 0.3;
  m_OmegaR = 8.4e-5;
  m_OmegaL = 0.7-m_OmegaR;
  m_Hubble = 0.7;          // H0 = m_Hubble*100 km/s/Mpc
  
}


void DIRBR_Redshift::InitializeIntegral(void)
{
  
  m_nsteps = ((int)1.0e5);
  m_qmin = 5.0e-4;
  m_grid_factor = -((double)m_nsteps)/log(m_qmin);
  m_dlnq = 1.0/m_grid_factor;
  
  double* farray = new double[m_nsteps];
  m_Farray = new double[m_nsteps];
  m_qarray = new double[m_nsteps];
  
  m_qarray[0] = 1.0;
  farray[0] = 0.0;
  m_Farray[0] = 0.0;
  for (int i = 1; i <m_nsteps; i++) {
    double y = exp(-double(i)/m_grid_factor);
    double factor = sqrt(1.0-y);
    farray[i] = 
      1.0/y*( (1.0+y-y*y/2.0)*log((1+factor)/(1-factor)) - (1.0+y)*factor ) ;
    m_qarray[i] = y;
  }
  
  for (int i = 1; i<m_nsteps; i++) {
    farray[i] = farray[i-1] + farray[i];
  }
  for (int i = 1; i<m_nsteps; i++) {
    m_Farray[i] = 2.0*m_qarray[i]*m_qarray[i]*m_qarray[i]*farray[i]*m_dlnq;
    //std::cout<<"q: "<<m_qarray[i]<<" F: "<<m_Farray[i]<<std::endl;
  }
  
  delete farray;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double DIRBR_Redshift::GetDIRBR( double lambda, double z)
/*!
  \note: IMPORTANT: to stay compatible with DIRBR.cpp class, and with the way
         the integration is done in PairProduction::PropagationEBL_redshift(),
	 we must return: Intensity(lambda,z)

  \return Intensity(lambda,z) Units: [nW m^-2 sr^-1]

*/
{
  
  /////////////////////////////////////////////////////////////////
  /// ERROR CHECKING: Make sure everything is in bounds.
  /////////////////////////////////////////////////////////////////
  if ((z<m_zmin) || (z>m_zmax) || (lambda<m_lambda_min) || (lambda>m_lambda_max) ) {
    std::cerr << "ERROR: ";
    if (z < m_zmin || z > m_zmax) {
      std::cerr<<"z outside redshift bounds.\n";
      std::cerr<<"z: "<<z<<std::endl;
      std::cerr<<"z_min: "<<m_zmin<<std::endl;
      std::cerr<<"z_max: "<<m_zmax<<std::endl;
    }
    else {
      std::cerr<<"lambda limit violated..."<<std::endl;
    }
    return 0.0;
  }

  double dirbr = SEDTable->LinInterpolate(lambda, z);
  return dirbr;
}


double DIRBR_Redshift::GetCMBR( double& lambda)
{
  
  double Tcmb = eV_K_B_d;  // T = 2.728 +/- 0.004 K
  Tcmb*=2.728;
  double e=1.2415/lambda;            // \mum -> eV
  
  return 5.042e+16*e*e*e*e/(exp(e/Tcmb)-1.);
  
}


double DIRBR_Redshift::OpticalDepth( double& E, double& RedShift )
/*! Calculates optical depth for a VHE photon traveling from 
    a source at redshift z. 

    \param E - detected by observer photon energy [TeV]

    \param z - redshift of a source [1]

    \return  optical depth  [1] 

*/
{
  
  //double OmegaM         = 0.3;
  //double OmegaL         = 0.7;

  double z = 0.0;
  double OpticalDepth = 0.0;

  double melec_ev = eV_MELEC_d;   // [eV]
  double hc = eVCM_HC_d;          // [eV cm]
  double kappa = hc/melec_ev/melec_ev;  // [eV^-1 cm]
  kappa *= E*1.0e16;                    // [mum] (eV -> TeV included...)

  double sigma_T = CGS_THOM_CS_d*1.0e-4;  // [m^2]
  double Ho = m_Hubble*100.0/kmToMpc_d;
  double TeV = 1.0e12;   // [eV]
  double eV_to_nJ = 1.60218e-10;
  double opt_depth_scale = 3.0/2.0*M_PI*sigma_T*(E*TeV)/melec_ev/melec_ev/Ho/eV_to_nJ;

// #ifdef TESTING
//   for (z=0.0; z < RedShift; z+=m_dz) {
    
//     double z1=(1.0+z);
//     double z2=z1*z1;
//     double z3=z1*z2;
//     double z4=z1*z3;

//     double Q = sqrt(OmegaM*z3+OmegaL);
    
//     double tau = 0.0;
//     for (int i = (m_nsteps-1); i>=0; i--) {
//       double lambda = kappa*z2*m_qarray[i];
//       if ( (z1*lambda) < m_lambda_min || (z1*lambda) > m_lambda_max) continue;
     
//       double dbr=GetDIRBR(lambda*z1,z);
//       //if(lambda > 100.) dbr+=GetCMBR(lambda);
//       tau += m_Farray[i]*dbr;
//     }

//     ///////////////////////////////////////////////////////
//     ////// WARNING: Missing cmbr in integration!!!!   /////
//     ///////////////////////////////////////////////////////
       
//     //tau *= (1.0/Q*m_dlnq);
//     tau *= z4/Q*m_dlnq;
    
//     if (m_dz < (RedShift-z)) OpticalDepth += tau*m_dz;
//     else OpticalDepth += tau*(RedShift-z);
    
//   }

// #else

  for (z=0.0; z < RedShift; z+=m_dz) {
    
    double z1=(1.0+z);
    double z2=(1.0+z)*z1;
    double z3=(1.0+z)*z2;

    double Q = sqrt(m_OmegaM*z3+m_OmegaL+(1.-m_OmegaM-m_OmegaL)*z2);
    
    double tau = 0.0;
    for (int i = (m_nsteps-1); i>=0; i--) {
      double lambda = kappa*z1*m_qarray[i];
      if (lambda < m_lambda_min || lambda > m_lambda_max) continue;
     
      double dbr=GetDIRBR(lambda,z);
      if(lambda > 100.) dbr+=GetCMBR(lambda);
      tau += m_Farray[i]*dbr;
    }
    tau *= z3/Q*m_dlnq;
    
    if (m_dz < (RedShift-z)) OpticalDepth += tau*m_dz;
    else OpticalDepth += tau*(RedShift-z);

  }

  OpticalDepth*=opt_depth_scale;

  return OpticalDepth;
  
}


double DIRBR_Redshift::Distance(double& z)
{


  return 0.0;
}



DIRBR_ERR DIRBR_Redshift::TestDIRBRlimits( void )
{
  // Nothing...Note sure if this NEEDS to be virtual?
  DIRBR_ERR error = SUCCESS;
  return error;
}


DIRBR_ERR DIRBR_Redshift::SetDIRBR(int Kwb, double * lambda_, double * I_ )
{
  // Nothing...Note sure if this NEEDS to be virtual?
  DIRBR_ERR error = SUCCESS;
  return error;
}
