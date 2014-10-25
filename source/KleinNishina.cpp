/*  KleinNishina.cpp
  KleinNishina class implementation file

  \author   Yusef Shafi            \n
            UCLA                   \n
	    yshafi@ucla.edu        \n

  \author   Timothy C. Arlen       \n
            UCLA                   \n
	    arlen@astro.ucla.edu   \n

  \author   Vladimir Vassiliev     \n
            UCLA                   \n
	    vvv@astro.ucla.edu     \n

  \author   Tom R. Weisgarber     \n
            UC Chicago            \n
	    trw@uchicago.edu      \n

  \date     08/17/2006

  \version  1.1

  \revision:
     03/11/2008 - Redshift Dependent Propagation Length
     08/26/2008 - Pass RelParticle directly, rather than its
	          components, so change functions for this.
     10/14/2008 - Changed the way delta_z is calculated for the electron in
                  function PropagationLengthBB_Redshift().

     6/5/2009   - IMPORTANT! Changed factor of 2 in calls to
                  IsotropicRadiation() and IsotropicSigma() as well as the
		  factor of 3/4 in IncompleteTotalCrossSection().

     3/16/2010  - Replacing all D0,D1, etc. with 0.0, 1.0, etc. since this
                  works now with the new version of qd/dd_real.h
  \note

*/


#include "KleinNishina.hpp"

using namespace PhysConst;

namespace IGCascade
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Overloaded class constructor
  /// \param _rng: rng
  KleinNishina::KleinNishina(TRandom3* _rng)
  {

    m_rng=_rng;

    m_DE = "1.0E-25";            // Relative computation precision of roots
    m_dx   = 5.0E-3;             // integration step from 0 to xU
    double xU   = 25.0;          // energy upper bound-physical units are xU*kT
    m_num_int = (int)(xU/m_dx);     // number of steps for Prop Length Integrals

    m_PI2d6=VEC3D_PI*VEC3D_PI/2.0/3.0;        // pi*pi/6 for PolyLog function

    m_egy_bb = "0.0";  // Initialize BB temp to zero   [eV]

    // For PropLength Integral:
    egy_min = "4.0e-6";  // [eV]
    egy_max = "10.0";   // [eV]
    m_num_bins_dec = "1000.0";
    m_egy_factor = pow(10.0, 1.0/m_num_bins_dec);
    VEC3D_T Energy = egy_min;
    while (Energy < egy_max) {
      m_lambda_vec.push_back(Double(1.0e4*eVCM_HC/Energy));
      Energy *= m_egy_factor;
    }

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  VEC3D_T KleinNishina::RelativisticKinematics(Vec4D& P_lep, Vec4D& P_p)
    /*! Computes relativistic kinematics of Compton scattering of massive
      particle and photon. Scattering angle is sampled utilizing
      full KleinNishina cross-section.

      \param  P_lep - Energy-Momentum 4-vector of massive particle,
                      in lab frame
              P_p - Energy-Momentum of massless photon in lab frame.

      \note: Returns relative difference between masses of a particle before
      and after interaction with photon (should always be zero)
    */
  {

    VEC3D_T m0 = sqrt(P_lep*P_lep);

    VEC3D_T gamma = (P_lep.r0)/m0;
    VEC3D_T gammabeta = sqrt(gamma*gamma - 1.0);
    Vec3D e = P_p.GetDirection();
    Vec3D n = P_lep.GetDirection();	 //!may return 0 if norm is 0

    VEC3D_T en = e*n;
    VEC3D_T root_e = sqrt(1.0 - en*en);   //!temporary parameter only used once
    Vec3D v = (e - en*n)/root_e; //!add random tiny EPS direction when root_e=0

    //Mandelstam invariant
    VEC3D_T s = (P_lep+P_p)*(P_lep+P_p);       //!s_=2E'/m
    VEC3D_T s_ = s/m0/m0 - 1.0;

    //Paramters a, b: perpendicular, parallel components of boosted e' vector
    VEC3D_T ab_divisor = gamma - gammabeta*en;
    VEC3D_T a = (gamma*en-gammabeta)/ab_divisor;
    VEC3D_T b = root_e/ab_divisor;

    // Sample y=(1-cos(theta))/2, where theta =
    // photon scattering angle in particle r.f.
    VEC3D_T y = Scattering(s_);
    // y=D0;

    //std::cout<<y<<" "<<std::endl;

    //compute sin_theta, cos_theta
    VEC3D_T sin_theta = 2.0 * sqrt(y * (1.0 - y));
    VEC3D_T cos_theta = 1.0 - 2.0 * y;

    //Generate random phi, calculate sin_phi, cos_phi
    VEC3D_T phi = (2.0 * VEC3D_PI)*((VEC3D_T) m_rng->Uniform());
    VEC3D_T sin_phi = sin(phi);
    VEC3D_T cos_phi = cos(phi);

    //Calculate E1'/E1, E1, e1
    VEC3D_T a_ = a*cos_theta + b*sin_theta*cos_phi;
    VEC3D_T b_ = b*cos_theta - a*sin_theta*cos_phi;
    VEC3D_T c = gamma + gammabeta*a_;
    VEC3D_T en1 = m0/2.0 * s_ * c/(1.0 + s_*y);
    Vec3D e1 = ((gamma*a_+gammabeta)*n + b_*v + sin_theta*sin_phi*(n^v))/c;

    //Update state of lepton, photon four momenta
    P_lep+=P_p;
    P_p = Vec4D(en1, en1*e1);
    P_lep-=P_p;

    s_ = P_lep*P_lep;

    //if ( s_ < D0 ) {
    if (s_ < 0.0) {
      std::cout<<"RelativisticKinematics: Loss of computation precision."<<
	std::endl;
      std::cout<<"RelativisticKinematics: Negative squared mass: "<<s_<<
	std::endl;
      exit(0);
    }
  //return D1 - sqrt(s_)/m0;
  return 1.0 - sqrt(s_)/m0;
    //return y;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  VEC3D_T KleinNishina::Scattering(VEC3D_T x)
    /* Samples Compton scattering angle of the photon in the particle's r.f.
       Returns y = (1-cos(theta))/2, such that no scattering (theta=0)
       corresponds to y=0.

       \param x = (s-m0^2)/m0^2 (Mandelstam invariant) = 2*E'/m0
                  (in particle's rest frame)

    */
  {

    VEC3D_T chi = (VEC3D_T) (m_rng->Uniform());

    // Klein Nishina total cross section:
    VEC3D_T KN_ = IncompleteTotalCrossSection(x,x);

    VEC3D_T LB = "0.0";			  // Left bound of the integral
    VEC3D_T RB = "1.0";			  // Right bound of the integral
    VEC3D_T y  = (LB + RB) / 2.0;
    VEC3D_T y_ = "0.0";

    VEC3D_T G;
    VEC3D_T z;

    while ( fabs(y - y_) > m_DE*fabs(y) ) {
      z = x*y;
      G = (IncompleteTotalCrossSection(z,x) - 3.0*y*y*(1.0 - y)/(1.0 + z))/KN_;
      y_ = y;

      if ( G < chi ) LB = y;
      else RB = y;

      y = (LB + RB) / 2.0;
    }
    return y;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  VEC3D_T KleinNishina::IncompleteTotalCrossSection(VEC3D_T z, VEC3D_T x)
    /*
      Klein_Nishina Function
      Computes the unitless Klein-Nishina total cross section
      from 0 to y, [y = (1-cos(theta))/2] given a value of z=xy
      and x (x = 2E/m=(s-m0^2)/m0^2). The total KN cross section is

      sigma (x) =  8/3*Pi*R_o^2 * IncompleteTotalCrossSection(x,x),

      where R_o is the classical electron radius.
    */
  {
    //if ( x <= "0.0" || z <= "0.0" ) {
    if ( x <= 0.0 || z <= 0.0 ) {
      std::cout<<"IncompleteTotalCrossSection: Non-positive argument."<<
	std::endl;
      exit(0);
      //return D0;
    }

    VEC3D_T arg = 1.0 + z;

    VEC3D_T Ln3;
    VEC3D_T Ln2;
    VEC3D_T Ln1;

    // Thomson regime
    if( z < 1.0/2.0 ) {

      int i = 2;
      Ln3  = "0.0";
      VEC3D_T dLn3 = "1.0";
      VEC3D_T zPower = -z * z;

      while ( fabs(dLn3) > m_DE*fabs(Ln3)) {
        i++;
	zPower = -zPower*z;
	dLn3=zPower/((VEC3D_T)((double)i));
	Ln3 += dLn3;
      }
      Ln2 = Ln3 - z*z/ 2.0;
      Ln1 = Ln2 + z;

    } else {

      // Klein-Nishina regime
      Ln1 = log(arg);
      Ln2 = Ln1-z;
      Ln3 = Ln2+z*z/2.0;
    }

    Ln1/=(x);
    Ln2/=(x*x);
    Ln3/=(x*x*x);

    return 3.0/4.0*(Ln1- 4.0*Ln2 - 8.0*Ln3 + z/x*(1.0 + arg)/2.0/arg/arg);

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  /////////////////////////////////////////////////////////////////////////////
  //NOTE: March 14, 2008-I commented out  last part of this function because //
  //       we can't convert VEC3D_T into a double or an int in any way that I//
  //       can find out. I need to figure out how to do this...              //
  //                                                                         //
  //UPDATE: Nov. 25, 2008-I reworked it using the header I wrote called      //
  //         convert.h. This converts VEC3D_T to double, which can then be   //
  //         converted to anything desired.                                  //
  /////////////////////////////////////////////////////////////////////////////

  bool KleinNishina::PropagationLengthStar(StructStar& star, Vec4D& P_lep,
			  Vec4D& R_e,Vec4D& P_p, Vec4D& R_p, VEC3D_T& pl)
    /*!
      Routine samples propagation length of a particle in
      a gas of stellar photons and samples parameters
      of a photon with which interaction takes place.

      \param Star - structure for an arbitrary star, containing:
             Radius - Radius of star				         [cm]
             Temperature - Star temperature			         [eV]
             Luminosity - Luminosity of star			         [eV/s]

      \param P_lep - 4-momentum of charged particle (charge e assumed)     [eV]
      \param R_e - 4-vector of position of charged particle                [cm]
      \param P_p - 4-momentum of stellar photon, interacting with particle [eV]
      \param R_p - 4-vector of position of charged particle                [cm]
      \param pl - the propagation length of charged particle	           [cm]

      \return  boolean value (true/false) for an interaction

      \note Some typical parameters for O, B stars:

                R/R_sun	     M/M_sun	L/L_sun		Temp (K)
      ------------------------------------------------------------
      O2	 16	      158    	200.E+4		5.40E+4
      05	 14	      58	80.E+4		4.60E+4
      B0         5.7	      16	1.6E+4		2.90E+4
      B5         3.7	      5.4	0.75E+4		1.52E+4
      A0         2.3	      2.6	0.063E+4	0.96E+4

      R_sun = 6.96E+10 cm;	M_sun = 1.988E+33 g;   L_sun = 2.389E+45 eV/sec

    */
  {

    // Define needed integration parameters.

    // Particle's parameters
    VEC3D_T m0 = sqrt(P_lep*P_lep);
    if ( m0 <= 0.0 ) {
      std::cout<<"PropagationLengthStar: Non-positive mass."<<std::endl;
      exit(0);
      //return;
    }

    VEC3D_T gamma = (P_lep.r0)/m0;
    if ( gamma <= 1.0 ) {
      std::cout<<"PropagationLengthStar: Invalid gamma."<<std::endl;
      exit(0);
      //return;
    }
    VEC3D_T gammabeta=sqrt( gamma*gamma-1.0 );

    VEC3D_T z = star.Temp/m0;			    // kT/mc^2 [no units]

    Vec3D n = P_lep.GetDirection();

    VEC3D_T rn = (VEC3D_T) (R_e.r*n);
    Vec3D Rho = R_e.r - rn*n;	                    // [cm]

    VEC3D_T magRho = Rho.Norm();
    VEC3D_T magRho2= magRho*magRho;
    VEC3D_T coseta_0 = rn/sqrt( magRho2 + rn*rn );
    VEC3D_T eta_0 = (acos(coseta_0));

    // Calculate factors of propagation length integral
    VEC3D_T IntPref = 4.0*VEC3D_PI * ( star.Radius*star.Radius/magRho);
    {
      VEC3D_T x=star.Temp/eVCM_HC;
      IntPref *= x*x*x;
    }
    IntPref *= CGS_THOM_CS*eV_MELEC*eV_MELEC/m0/m0;


    //-------------------------------------------------------
    // Generate Random number between 0 & 1
    VEC3D_T chi = (VEC3D_T) (m_rng->Uniform());
    VEC3D_T LnChi = - log(chi);
    VEC3D_T eta = eta_0;
    VEC3D_T dEta = (VEC3D_T)(VEC3D_PI/10000.);	// eta integration step

    VEC3D_T * IS = new VEC3D_T [m_num_int];
    int j;

    VEC3D_T tau = "0.0";
    VEC3D_T dtau;
    VEC3D_T Sum_IS = "0.0";
    while ( ( tau < LnChi ) && ( eta > 0.0 ) ) {
      VEC3D_T coseta=cos(eta);
      VEC3D_T S = "0.0";
      VEC3D_T x = m_DE;
      VEC3D_T u;

      for (j = 0; j < m_num_int; j++) {
	u = z*x*(gamma - gammabeta * coseta);
	S += (x*x / ( exp(x) - 1.0))* IncompleteTotalCrossSection(u, u);
	IS[j]=S;
	x+= m_dx;
      }
      Sum_IS = S;

      S *= fabs( gammabeta - gamma*coseta )/gammabeta;
      S *= IntPref*m_dx;

      dtau = S*dEta;
      tau += dtau;

      eta -= dEta;
    }

    if ( eta < 0.0 ) {  // No interaction, particle leaves stellar photon field
      delete IS;
      return false;
    }
    // interpolation scheme
    eta -=dEta*(tau-LnChi)/dtau;
    pl = magRho/tan(eta) - rn;

    std::cout<<"\n eta = "<<eta<<std::endl<< "eta_0 = " << eta_0 << std::endl
	      << " n = " << n << std::endl << std::endl;


    // Determine energy of photon that charged particle interacts with
    VEC3D_T chi_2 = (VEC3D_T) (m_rng->Uniform());
    //    VEC3D_T PhotEnerInt = D0;
    //--------------------------------------------

    VEC3D_T IntegralRatio = "2.0";
    VEC3D_T RB = (VEC3D_T) m_num_int;
    VEC3D_T LB = "0.0";
    VEC3D_T Steps = "0.0";
    int IntSteps = 0;
    VEC3D_T CompParam = 0.05;		// Need to obtain more accuracy!

    j=0;
    while ( fabs(chi_2 - IntegralRatio) > CompParam ) {
      Steps = (RB - LB)/2.0;
      // Testing:
      if (j == 0 || j == 1 || j == 2 || j == 3){
	std::cout << "j = " << j << "Steps = " << IntSteps << std::endl;
      }

      IntSteps = (int) Double(Steps);
      IntegralRatio = IS[IntSteps]/IS[m_num_int-1];
      //      IntegralRatio = IS[IntSteps]/IS[m_num_int-1];
      //      	std::cout<<"IntegralRatio = "<<IntegralRatio<<std::endl;

      if ( IntegralRatio < chi_2 ) LB = Steps;
      else RB = Steps;
      j++;
      std::cout<<"fabs(chi_2 - IntegralRatio) = "<<fabs(chi_2 - IntegralRatio)
	       <<std::endl;

    }
    //---------------------------------------------


    /*
      for( j = m_num_int - 100; j < m_num_int - 1; j++)
      {
      std::cout << "IS[j]/IS[m_num_int-1] = " << IS[j]/IS[m_num_int-1] << std::endl;
      }
    */

    /*
      j = 0;
      while ( ( PhotEnerInt < chi_2 ) && ( j < m_num_int ) ) {
      PhotEnerInt = IS[j]/Sum_IS;
      j++;
      }
    */


    VEC3D_T x = ( Steps - 1.0 )*m_dx;  // holds the energy of incoming photon
    std::cout << std::endl << std::endl << "chi_2 = " << chi_2 << std::endl <<
      "IntegralRatio = " << IntegralRatio <<std::endl<<"Energy of Photon = "<<
      x << " kT." << std::endl;

    delete IS;

    // Update Position of particle and photon
    R_e.r = (R_e.r + n*pl)/star.Radius;
    R_p.r = R_e.r;

    // Update 4-momentum of photon
    P_p.r0 = x * star.Temp;
    P_p.r = P_p.r0 * R_p.r;

    //	std::cout<<"Photon final position = "<<R_p.r<<std::endl<<
    // "Particle Final Position = "<<R_e.r << std::endl;
    //-------------------------------------------------------

    return true;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  VEC3D_T KleinNishina::PropagationLengthBB(VEC3D_T  Ebb, Vec4D& P_lep,
					    Vec4D& P_p)

  /*!
    Routine samples propagation length of a particle in
    the gas of photons with BB spectrum and samples 4 momentum
    of the photon with which interaction takes place.

    \param Ebb   - kT of BB spectrum                                       [eV]
    \param P_lep - 4-momentum of charged particle (charge e assumed),
                   is NOT modified in this function.                     [eV]
    \param P_p   - 4-vector of BB photon, interacting with particle,
                   updated at end of function                            [eV]

    \return propagation length                                           [cm]

  */

  {
    // Particle's parameters
    VEC3D_T m0 = sqrt(P_lep*P_lep);
    if ( m0 <= 0.0 ) {
      std::cout<<"PropagationLengthBB: Non-positive mass."<<std::endl;
      exit(0);
      //return;
    }
    //std::cout<<"m0 = "<<m0<<std::endl;

    VEC3D_T gamma = (P_lep.r0)/m0;
    if ( gamma <= 1.0 ) {
      std::cout<<"PropagationLengthBB: Invalid gamma."<<std::endl;
      exit(0);
      //return;
    }

    // Compute BB parameters
    if( fabs(1.0-m_egy_bb/Ebb) > m_DE ) {
      VEC3D_T q=Ebb/eVCM_HC;
      q=4.0*VEC3D_PI*q*q*q;
      m_num_dens_bb=q*2.0*(2.0*1.2020569032); // Photon density [cm^-3]
      m_egy_dens_bb=q*m_PI2d6*m_PI2d6*8.0*3.0/5.0*Ebb;  // Photon energy density [eV cm^-3]
      m_egy_bb=Ebb;

      //VEC3D_T mfp_=eV_MELEC/CGS_THOM_CS/m_egy_dens_bb*3.0/4.0;  // units??
      //std::cout<<"mfp_: "<<mfp_<<std::endl;
      //std::cout<<"m_num_dens_bb: "<<m_num_dens_bb<<std::endl;
      //std::cout<<"m_egy_dens_bb: "<<m_egy_dens_bb<<std::endl;
    }

    VEC3D_T z=Ebb/m0;

    // Compute BB spectrum integral
    VEC3D_T * IS = new VEC3D_T [m_num_int];
    int j;
    VEC3D_T x = m_DE;
    VEC3D_T S = "0.0";
    for (j=0; j<m_num_int; j++) {
      IS[j]=x*x / (exp(x)-1.0) * IsotropicSigma(2.0*z*x,gamma);
      S+=IS[j];
      x+=m_dx;
    }
    S+=-(IS[0]+IS[m_num_int-1])/2.0;
    S*=m_dx;

    // Mean free path
    VEC3D_T mfp;
    {
      VEC3D_T ST = CGS_THOM_CS*eV_MELEC*eV_MELEC/m0/m0;
      VEC3D_T q=Ebb/eVCM_HC;
      q=4.0*VEC3D_PI*q*q*q;
      mfp=1.0/(q*ST*S);
      //std::cout<<"mfp BB = "<<mfp<<std::endl;
    }

    // Sample Photon's energy
    VEC3D_T chi = (VEC3D_T) (m_rng->Uniform());
    VEC3D_T S_ = "0.0";
    j=0;
    x = "0.0";
    while (chi > S_/S && j < m_num_int-1) {
      S_+=m_dx*(IS[j]+IS[j+1])/2.0;
      x+=m_dx;
      j++;
    }
    x-=2.0*(S_- S*chi)/(IS[j]+IS[j-1]);

    delete[] IS;

    // Photon's energy
    VEC3D_T E_p = Ebb * x;

    // Sample cos(theta) and phi
    VEC3D_T cs = IsotropicRadiation(2.0*z, gamma);
    VEC3D_T sn = sqrt(1.0 - cs*cs);
    VEC3D_T phi = (2.0 * VEC3D_PI)*((VEC3D_T) m_rng->Uniform());
    VEC3D_T sin_phi = sin(phi);
    VEC3D_T cos_phi = cos(phi);
    // std::cout<<"cos(theta): "<<cs<<std::endl;

    // Define basis
    Vec3D n = P_lep.GetDirection();
    Vec3D e1 = Vec3D(n.y, -n.x, 0.0);
    VEC3D_T nrm = e1.Norm();
    if (nrm == 0.0) e1=Vec3D(1.0,0.0,0.0);
    else e1/=nrm;
    Vec3D e2 = n^e1;

    // Determine photon's Momentum Vec4
    Vec3D e = cs*n + sn*sin_phi*e2 + sn*cos_phi*e1;
    P_p=Vec4D(E_p,E_p*e);
    // std::cout<<"mass2: "<<P_p*P_p<<std::endl;

    VEC3D_T chi_ran = (VEC3D_T) (m_rng->Uniform());
    //VEC3D_T chi_ran = 1.E-8;
    return -mfp*log(chi_ran);

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  VEC3D_T KleinNishina::PropagationLengthCMBandIR(DIRBRBase* ebl_model,
			  VEC3D_T Ebb, Vec4D& P_lep, Vec4D& P_p)
    /*!
      Same as PropagationLengthBB() except it adds the IR Energy Density of
      the ebl spectrum to sample the Propagation Length and energy of
      upscattered photon.

      \param ebl_model - ebl model to use.
      \param Ebb - kT of BB spectrum                                       [eV]
      \param P_lep - 4-momentum of charged particle (charge e assumed),
                     is NOT modified in this function.                     [eV]
      \param P_p   - 4-vector of BB photon, interacting with particle,
                     updated at end of function                            [eV]

      \return propagation length                                           [cm]

    */
  {

    // Particle's parameters
    VEC3D_T m0 = sqrt(P_lep*P_lep);
    if ( m0 <= 0.0 ) {
      std::cout<<"PropagationLengthBB: Non-positive mass."<<std::endl;
      exit(0);
      //return;
    }
    //std::cout<<"m0 = "<<m0<<std::endl;

    VEC3D_T gamma = (P_lep.r0)/m0;
    if ( gamma <= 1.0 ) {
      std::cout<<"PropagationLengthBB: Invalid gamma."<<std::endl;
      exit(0);
      //return;
    }

    // Compute Photon field spectrum integral
    VEC3D_T epsilon = egy_min;
    VEC3D_T delta_epsilon = epsilon;
    std::vector<VEC3D_T> IS;
    std::vector<VEC3D_T> Delta_epsilon;
    VEC3D_T S = "0.0";
    unsigned Nstep = 0;
    VEC3D_T CMB_max = 25.0*Ebb;
    while(epsilon < egy_max) {

      VEC3D_T temp_dirbr = ebl_model->GetDIRBR(m_lambda_vec[Nstep],0.0);
      if(epsilon < CMB_max) temp_dirbr+=ebl_model->GetCMBR(m_lambda_vec[Nstep]);

      VEC3D_T temp_IS = temp_dirbr*IsotropicSigma(2.0*epsilon/m0,gamma)
	/epsilon/epsilon;

      IS.push_back(temp_IS);
      S+=temp_IS*delta_epsilon;
      Delta_epsilon.push_back(delta_epsilon);

      VEC3D_T epsilon_new = epsilon*m_egy_factor;
      delta_epsilon = epsilon_new - epsilon;
      epsilon = epsilon_new;
      Nstep++;
    }

    // Mean free path

    VEC3D_T mfp;
    {
      VEC3D_T ST = CGS_THOM_CS*eV_MELEC*eV_MELEC/m0/m0;
      mfp=1.0/(ST*S*2.0*VEC3D_PI/CGS_C*nJ_TO_eV*1e-4);
      //std::cout<<"mfp TA = "<<mfp<<std::endl;
    }

    // Sample Photon's energy
    VEC3D_T chi = (VEC3D_T) (m_rng->Uniform());
    VEC3D_T S_ = "0.0";
    unsigned j=0;
    epsilon = egy_min;
    while (chi > S_/S && j < Nstep) {
      S_+=IS[j]*Delta_epsilon[j];
      j++;
      epsilon *= pow(10,1.0/m_num_bins_dec);
    }

    // Photon's energy
    VEC3D_T E_p = epsilon;
    //std::cout<<"Energy of Bg Photon: "<<E_p<<std::endl;

    // Sample cos(theta) and phi
    VEC3D_T cs = IsotropicRadiation(2.0*Ebb/m0, gamma);
    VEC3D_T sn = sqrt(1.0 - cs*cs);
    VEC3D_T phi = (2.0 * VEC3D_PI)*((VEC3D_T) m_rng->Uniform());
    VEC3D_T sin_phi = sin(phi);
    VEC3D_T cos_phi = cos(phi);
    // std::cout<<"cos(theta): "<<cs<<std::endl;

    // Define basis
    Vec3D n = P_lep.GetDirection();
    Vec3D e1 = Vec3D(n.y, -n.x, 0.0);
    VEC3D_T nrm = e1.Norm();
    if (nrm == "0.0") e1=Vec3D(1.0,0.0,0.0);
    else e1/=nrm;
    Vec3D e2 = n^e1;

    // Determine photon's Momentum Vec4
    Vec3D e = cs*n + sn*sin_phi*e2 + sn*cos_phi*e1;
    P_p=Vec4D(E_p,E_p*e);

    VEC3D_T chi_ran = (VEC3D_T) (m_rng->Uniform());
    //VEC3D_T chi_ran = 1.E-8;
    return -mfp*log(chi_ran);

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  VEC3D_T KleinNishina::IsotropicSigma(VEC3D_T z, VEC3D_T gamma)
    /*
      This routine computes photon direction averaged cross-section
      weighted by the relative velocity of particle and photon

      <sigma * |v-c*cos(theta)|> / v

      \param  z= E/mc^2  energy of the photon in lab r.f. divided by
      particle's mass
      \param  gamma      particle's gamma

      \note To get correct physical units the result must
      be multiplied by  Thomson cross-section 8/3*Pi*R_o^2
    */
  {
    if ( gamma <= 1.0 || z <= 0.0 ) {
      std::cout<<"IsotropicSigma:  Ivalid argument(s)."<<std::endl;
      exit(0);
      //return;
    }

    if ( z < m_DE ) z=m_DE;   // prevents loss of computational accuracy

    VEC3D_T gammabeta = sqrt(gamma*gamma - 1.0);
    VEC3D_T b=1.0/gamma;

    if( fabs(1.0-m_IRd_b/b) > m_DE || fabs(1.0-m_IRd_z/z) > m_DE ) {
      VEC3D_T a=gamma+gammabeta;
      VEC3D_T c=1.0/a;
      IsotropicFactors(m_IRd_f1c, m_IRd_f2c, z*c);
      IsotropicFactors(m_IRd_f1b, m_IRd_f2b, z*b);
      IsotropicFactors(m_IRd_f1a, m_IRd_f2a, z*a);
      m_IRd_U=-b*(m_IRd_f1a + m_IRd_f1c - 2.0*m_IRd_f1b)/z
	+ (m_IRd_f2a+m_IRd_f2c-2.0*m_IRd_f2b)/z/z;
      m_IRd_b=b;
      m_IRd_z=z;
    }

    return gamma*m_IRd_U/gammabeta/gammabeta/gammabeta;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  VEC3D_T KleinNishina::IsotropicRadiation(VEC3D_T  z, VEC3D_T gamma)
    /*
      Routine samples impact angle of the photon that interacts with
      the particle. In the lab reference frame the photon field is
      isotropic. This distribution is changed by the relative
      velocity factor as well as dependence of cross-section
      on photon's energy (KN regime)

      sigma * |v-c*cos(theta)| / v

      \param  z= E/mc^2  energy of the photon in lab r.f. divided by
      particle's mass
      \param  gamma      particle's gamma


    */
  {
    if ( gamma <= 1.0 || z <= 0.0 ) {
      std::cout<<"IsotropicRadiation: Invalid argument."<<std::endl;
      exit(0);
      //return;
    }

    VEC3D_T b=1.0/gamma;
    //if( fabs(1.0-m_IRd_b/b) > m_DE || fabs(1.0-m_IRd_z/z) > m_DE ) {
      //VEC3D_T S=IsotropicSigma(z, gamma);
    //}

    VEC3D_T chi = (VEC3D_T) (m_rng->Uniform());

    VEC3D_T gammabeta = sqrt(gamma*gamma - 1.0);
    VEC3D_T RB = gamma+gammabeta;		 // Right bound of the integral
    VEC3D_T LB = 1.0/RB;			  // Left bound of the integral
    VEC3D_T y  = (LB + RB)/2.0;
    VEC3D_T y_ = LB;

    VEC3D_T UI;
    VEC3D_T f1;
    VEC3D_T f2;

    while ( fabs(y - y_) > m_DE*fabs(y) ) {

      IsotropicFactors(f1, f2, z*y);
      if (y <= b) UI=m_IRd_b*(f1-m_IRd_f1c)/z-(f2-m_IRd_f2c)/z/z;
      else {
	UI = -m_IRd_b*(f1+m_IRd_f1c-2.0*m_IRd_f1b)/z+
	  (f2+m_IRd_f2c-2.0*m_IRd_f2b)/z/z;
      }
      y_ = y;
      if (UI/m_IRd_U<chi) LB = y;
      else RB = y;
      y = (LB + RB)/2.0;
    }

    return (gamma-y)/gammabeta;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void KleinNishina::IsotropicFactors(VEC3D_T & f1, VEC3D_T & f2, VEC3D_T x)
    /*
      Klein_Nishina Isotropic factors
      Factors f1 and f2 are required to sample impact angle of photon
      and particle. This routine assumes isotropic distribution of
      photons in the lab frame.
    */
  {
    if ( x < 0.0 ) {
      std::cout<<"IsotropicFactors: Negative argument."<<std::endl;
      exit(0);
      //f1=D0;
      //f2=D0;
      //return;
    }

    VEC3D_T arg = 1.0 + x;
    VEC3D_T Ln1;
    VEC3D_T Ln2;
    VEC3D_T Ln3;

    if( x < 1.0/2.0 ) {  // small x regime
      Ln3  = "0.0";

      int i = 2;
      VEC3D_T dLn3 = 1.0;
      VEC3D_T xPower = -x * x;
      while ( fabs(dLn3) > m_DE*fabs(Ln3)) {
	i++;
	xPower = -xPower*x;
	dLn3=xPower/((VEC3D_T)((double)i));
	Ln3 += dLn3;
      }
      Ln2 = Ln3 - x*x / 2.0;
      Ln1 = Ln2 + x;

    } else {           // large x regime
      Ln1 = log(arg);
      Ln2 = Ln1-x;
      Ln3 = Ln2+x*x/2.0;
    }

    VEC3D_T PLg = PolyLog1(x);

    if (x>0.0) {
      Ln2/=x;
      Ln3/=(x*x);
    }
    f1=3.0/4.0 * (x/arg/2.0 + Ln1/2.0 + 4.0*Ln2 + 4.0*Ln3 + PLg);
    f2=3.0/4.0*(arg*Ln1 - x*(1.0+arg)/arg/2.0 + 8.0*Ln1 + 8.0*Ln2 - 4.0* PLg);

    return;
  }



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ///////////////////////////////////////////////////////////////////
  // ABSTRACT MATHEMATICAL FUNCTIONS ////////////////////////////////
  // \note: should be moved /////////////////////////////////////////
  // \note: not directly relevant to the KleinNishina class /////////
  ///////////////////////////////////////////////////////////////////

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  VEC3D_T KleinNishina::PolyLog1(VEC3D_T x)
    /*
      Dilogarithm related function. In terms of traditional notation for
      dilogarithm, Li_2(z)=z+z^2/4+z^3/9+...+z^n/n^2, it is defined as

      PolyLog1(x)=-Li_2(-x)

      In terms of integral representation

      PolyLog1(x)=Int_1^x  ln(1+x)m_dx/x+pi^2/12

      The function is defined for x from 0 to infinity and PolyLog1(1)=pi^2/12

    */
  {
    if ( x < 0.0 ) {
      std::cout<<"PolyLog1: Negative argument."<<std::endl;
      exit(0);
      //VEC3D_T PL=D0;
      //return PL;
    }

    if ( x > 2.0 ) {              // large x approximation  (x>2)

      VEC3D_T Ln=log(x);
      VEC3D_T xPower = "1.0";
      VEC3D_T PL  = Ln*Ln/2.0 + m_PI2d6;
      VEC3D_T dPL = m_PI2d6;
      VEC3D_T i = "0.0";

      while (fabs(dPL) > m_DE*m_PI2d6 ) {
	i+=1.0;
	xPower/=(-x);
	dPL=xPower/i/i;
	PL+=dPL;
      }
      //std::cout<<" 1 "<<PL<<std::endl;
      return PL;

    } else if ( x < ((VEC3D_T) 0.618)) {       // small x approximation

      VEC3D_T xPower = "1.0";
      xPower = -xPower;
      VEC3D_T PL  = "0.0";
      VEC3D_T dPL = m_PI2d6;
      VEC3D_T i = "0.0";

      while (fabs(dPL) > m_DE*m_PI2d6 ) {
	i+=1.0;
	xPower*=(-x);
	dPL=xPower/i/i;
	PL+=dPL;
      }
      //std::cout<<" 2 "<<PL<<std::endl;
      return PL;
    }

    VEC3D_T xPower = "1.0";                      // x~1 regime
    xPower = -xPower;

    VEC3D_T Ln=log( 1.0 + x );
    VEC3D_T PL  = log (x)*Ln - Ln*Ln/2.0 +  m_PI2d6;
    VEC3D_T dPL = m_PI2d6;
    VEC3D_T i = "0.0";
    VEC3D_T y = 1.0 / (1.0 + x);

    while (fabs(dPL) > m_DE*m_PI2d6 ) {
      i+=1.0;
      xPower*=y;
      dPL=xPower/i/i;
      PL+=dPL;
    }
    //std::cout<<" 3 "<<PL<<std::endl;
    return PL;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}


