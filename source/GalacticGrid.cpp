/*!
-------------------------------------------------------------------------------
  \file     GalacticGrid.cpp
            Propagates charged leptons and photons
	    through cosmological magnetic fields.

  \author   Timothy C. Arlen                     \n
            Department of Physics and Astronomy  \n
            UCLA                                 \n
            arlen@astro.ucla.edu                 \n

  \author   Yusef Shafi                          \n
            UCLA	                         \n
	    yshafi@ucla.edu	                 \n

  \author   Stephen Fegan                        \n
            UCLA                                 \n
	    sfegan@astro.ucla.edu                \n

  \date     03/14/2008

  \version  1.1

  \revision 04/24/2008 - Added Magnetic Field Propagation in Expanding universe
                         in PropagateMagFieldExpansion() function.
	    02/06/2009 - Completely changed the way I determined if
	                 propagation would take electron to a new cell.

  \note     1) IMPORTANT:
            This function only works properly for
	      - cell_size >= 0.1 Mpc, and for
	      - B >= 10E-10
	      - Electron Energy > 0.1 TeV
	    due to issues in correctly computing cell-crossing for electron IC
	    propagation steps in the cmb photon field.
-------------------------------------------------------------------------------
*/

#include "GalacticGrid.hpp"

using namespace PhysConst;

namespace IGCascade
{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Overloaded class constructor (1)
  /// \param _MagneticField: m_MagneticField
  /// \param _rng: rng
  MagneticGrid::MagneticGrid(TRandom3 * _rng, VEC3D_T B_mag,
			     std::string s_cell_size, std::string filename)
  {

    // public member:
    m_DE = "1.0E-25";

    //Field _MagneticField;
    //m_MagneticField = _MagneticField;

    m_rng = _rng;

    // cellsize defined in units of Mpc but converted to cm
    std::istringstream(s_cell_size) >> m_cellsize; // [Mpc]
    //MpcToCm = 3.086E+24;                         // [cm/Mpc]
    m_cellsize = m_cellsize*Double(MPC_TO_CM);     // [cm]
    m_bmag = B_mag;                                // [gauss]
    m_sfilename = filename;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void MagneticGrid::
  PropagateBFieldRedshift(RelParticle& Photon,
			  RelParticle*& Lepton, Vec3D& n_eo, VEC3D_T& PL,
			  VEC3D_T& delta_z, const bool LOCK)
  /* Designed to work with the KleinNishina class. A propagation
     length has been determined and the IC scattered photon and
     lepton are computed, then the lepton and photon are translated
     and rotated due to the presence of the cosmological magnetic field.

     \param:

       Photon - The photon which will be rotated and translated

       Lepton - The lepton which will be rotated and
	        translated

       n_eo - Lepton direction before IC scat with photon.

       PL - Propagation length of lepton                     [cm]

       delta_z - change in redshift of the lepton over the
                 distance of the propagation length.

       LOCK - if 'true', use CheckMagneticField_Lock; default value is 'false'.

       IMPORTANT: assumes delta_z is always a positive quantity.

      \note: returns void, but updates photon/lepton 4-position and
             momentum directly
    */
  {

    VEC3D_T PL_original = PL;
    VEC3D_T MIN_Precision = 1.E-12*PL_original;
    VEC3D_T PL_remain = PL;
    while (PL_remain > MIN_Precision) {

      PL_remain = PL;

      ICoord c_cur = CheckCurrentCell(Lepton->m_r4.r);

      Vec3D e_b(0.,0.,0.);
      if (LOCK) e_b = CheckMagneticField_Lock(c_cur);
      else e_b = CheckMagneticField(c_cur);

      Vec3D r_new = UpdatePosition(Lepton, PL, n_eo, e_b);
      ICoord c_new = CheckCurrentCell(r_new);

      if (c_new != c_cur) {  // Cell crossed
	//std::cerr<<"Cell Crossing"<<std::endl;
	// Find Length to new cell:
	VEC3D_T PL_left = 0;
	VEC3D_T PL_right = PL;
	while( (PL_right - PL_left) > MIN_Precision ) {
	  VEC3D_T PL_test = 0.5*(PL_left+PL_right);
	  r_new = UpdatePosition(Lepton, PL_test, n_eo, e_b);
	  ICoord c_test = CheckCurrentCell(r_new);
	  if(c_test == c_cur) PL_left = PL_test;
	  else PL_right = PL_test;
	}
	//VEC3D_T PL_fraction = PL_right/PL;
	PL = PL_right;
      }

      VEC3D_T delta_z  = Delta_z(PL, Lepton);
      VEC3D_T delta_time = 0.0;
      VEC3D_T delta_zs = Delta_zs(PL,Lepton,r_new,delta_time);

      if ( (Lepton->m_z_s - delta_zs) >= m_DE ) { // z_s = 0 surface not crossed

	Propagation(Photon,Lepton,PL,r_new,delta_time,n_eo,
		    delta_z,delta_zs,e_b);
	PL_remain -= PL;

      }
      else {  // z_s = 0 surface crossed

	/////////////////////////////////////////////
	// compute PL to z_s = 0 surface:          //
	/////////////////////////////////////////////
	VEC3D_T PL_right = PL;
	VEC3D_T PL_left = "0.0";
	Vec3D r_e = Lepton->m_r4.r;
	// Compute PL to z_s = 0.0 surface.
	while ( fabs(Lepton->m_z_s - delta_zs) > m_DE ) // About 40 iterations.
	  {
	    PL = (PL_right + PL_left)/2.0;
	    r_new = UpdatePosition(Lepton, PL, n_eo, e_b);
	    delta_zs = Delta_zs(PL,Lepton,r_new,delta_time);

	    if ( (Lepton->m_z_s - delta_zs) > 0.0 ) PL_left = PL;
	    else PL_right = PL;
	  }
	/////////////////////////////////////////////

	delta_z = Delta_z(PL, Lepton);

	Propagation(Photon,Lepton,PL,r_new,delta_time,n_eo,
		    delta_z,delta_zs, e_b);
	PL_remain = "0.0";

      }

    } // end while


  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  Vec3D MagneticGrid::UpdatePosition(RelParticle*& Lepton, VEC3D_T& PL,
				     Vec3D& n_eo, Vec3D& e_b)
    /*!
      Computes the new position of Lepton if propagated a distance PL.
        IMPORTANT: does not change any of Lepton's values or PL.

      \param PL     - prop length of lepton [cm]
      \param Lepton - charged lepton which is to be propagated

      \return r_new - new position of Lepton
    */

  {

    VEC3D_T ze   = Lepton->m_z;
    //VEC3D_T m_eV = sqrt(Lepton->m_p4*Lepton->m_p4);
    //std::cerr<<"me_eV (Update): "<<m_eV<<std::endl;
    VEC3D_T p_o  = Lepton->m_p4.r.Norm()/(1.0+ze);
    //VEC3D_T rho  = m_eV*m_eV/p_o/p_o;              // [unitless]

    VEC3D_T Kappa_o = CGS_C*1.E-8;
    Kappa_o *= Lepton->m_q*m_bmag/p_o;

    VEC3D_T KappaL = Kappa_o*PL;

    Vec3D A = n_eo*PL;
    Vec3D DeltaR = A;
    // Play with m_bmag values to see what best precision is:
    if (fabs(KappaL) > 0.0)
      {
	Vec3D B = (n_eo - (n_eo*e_b)*e_b)*(sin(KappaL) - KappaL)/Kappa_o;
	Vec3D C = (n_eo^e_b)*(1.0 - cos(KappaL))/Kappa_o;
	DeltaR += (B + C);
      }

    Vec3D r_new = Lepton->m_r4.r + DeltaR;

    return r_new;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  VEC3D_T MagneticGrid::Delta_z(VEC3D_T& PL, RelParticle*& Lepton)
    /*!
      Calculates change in redshift for a lepton propagating a distance PL
        IMPORTANT: does not change any of Lepton's values.

      \param PL     - prop length calculated by PropagationLengthBB [cm]
      \param Lepton - charged lepton which is created by IC process

      \return delta_z - change in redshift of lepton

      **NOTE** Changed: E -> E_i *( 1+z_{i+1} )/( 1+z_i) in redshift expansion
    */
  {

    VEC3D_T m_eV = sqrt(Lepton->m_p4*Lepton->m_p4);
    VEC3D_T ze = Lepton->m_z;
    VEC3D_T p_o = Lepton->m_p4.r.Norm()/(1.0 + ze);
    VEC3D_T rho = m_eV*m_eV/p_o/p_o;

    VEC3D_T z1 = (1.0+ze);
    VEC3D_T z2 = z1*z1;
    VEC3D_T z3 = z1*z2;
    VEC3D_T z4 = z1*z3;
    VEC3D_T Q=sqrt(OMEGA_R*z4+OMEGA_M*z3+OMEGA_L+(1.0 - OMEGA_0)*z2);

    VEC3D_T F = 1.0/Q;
    VEC3D_T F2 = F*F;
    VEC3D_T F3 = F*F2;
    VEC3D_T G = 2.0*OMEGA_R*z3 + 3.0/2.0*OMEGA_M*z2 + (1.0-OMEGA_0)*z1;
    VEC3D_T H = 6.0*OMEGA_R*z2 + 3.0*OMEGA_M*z1 + (1.0-OMEGA_0);
    //VEC3D_T J  = 12.0*OMEGA_R*z1 + 3.0*OMEGA_M;
    VEC3D_T F_pr1= -F3*G;
    VEC3D_T F_pr2= -3.0*F2*F_pr1*G - F3*H;
    //VEC3D_T F_pr3= -3.0*F*(2.0*F_pr1*F_pr1 + F*F_pr2)*G - 6.0*F2*F_pr1*H -F3*J;

    VEC3D_T P = (1.0 - rho/z2/2.0);
    VEC3D_T P_pr1 = rho/z3;
    VEC3D_T P_pr2 = -3.0*rho/z4;
    //VEC3D_T P_pr3 = 12.0*rho/(z2*z3);

    VEC3D_T a = P*F;
    VEC3D_T a2 = a*a;
    VEC3D_T a3 = a*a2;
    //VEC3D_T a4 = a*a3;
    VEC3D_T b = P_pr1*F + P*F_pr1;
    VEC3D_T b2 = b*b;
    //VEC3D_T b3 = b*b2;
    VEC3D_T c = P_pr2*F + 2.0*P_pr1*F_pr1 + P*F_pr2;
    //VEC3D_T d = P_pr3*F + 3.0*P_pr2*F_pr1 + 3.0*P_pr1*F_pr2 + P*F_pr3;

    VEC3D_T A = 1.0/a;
    VEC3D_T B = b/2.0/a3;
    VEC3D_T C = (b2/2.0/a2 - c/6.0/a)/a3;
    //VEC3D_T D = (d/24.0 - 5.0/12.0*b*c/a + 5.0/8.0*b3/a2)/a3/a2;

    VEC3D_T L_prop = PL/CGS_HUBRAD;
    VEC3D_T L_prop2 = L_prop*L_prop;
    VEC3D_T L_prop3 = L_prop*L_prop2;
    //VEC3D_T L_prop4 = L_prop*L_prop3;

    VEC3D_T delta_z = A*L_prop + B*L_prop2 + C*L_prop3;// + D*L_prop4;

    return delta_z;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  VEC3D_T MagneticGrid::Delta_zs(VEC3D_T& PL, RelParticle*&
				 Lepton, Vec3D& r_new, VEC3D_T& delta_time)
    /*!
      Calculates delta_zs of a radial PHOTON AND time_delay given prop length
        PL.

      IMPORTANT:
        - does not change any of Lepton's values
	- Modifies delta_time

      \param PL     - Length to be propagated by Electron.
      \param Lepton - charged lepton which is created by IC process.
      \param r_new  - updated position of Lepton.

      \return delta_z_s - change in redshift of lepton

      \note  -  IMPORTANT: I couldn't exactly verify that the 4th order
                correction is correct or does anything, for some reason. So
		we just go up to third order on time computation.

    */
  {

    VEC3D_T m_eV = sqrt(Lepton->m_p4*Lepton->m_p4);
    //std::cerr<<"m_eV (Delta_zs): "<<m_eV<<std::endl;
    VEC3D_T ze = Lepton->m_z;
    VEC3D_T p_o = Lepton->m_p4.r.Norm()/(1.0 + ze);
    VEC3D_T rho = m_eV*m_eV/p_o/p_o;

    //------------------------ BEGIN z_s Computation ------------------------//
    VEC3D_T z = Lepton->m_z_s;
    VEC3D_T z1 = (1.+z);
    VEC3D_T z2 = z1*z1;
    VEC3D_T z3 = z2*z1;
    VEC3D_T z4 = z3*z1;
    VEC3D_T Q  = sqrt(OMEGA_R*z4+OMEGA_M*z3+OMEGA_L+(1.0-OMEGA_0)*z2);

    VEC3D_T a = 1.0/Q;
    VEC3D_T a2 = a*a;
    VEC3D_T a3 = a*a2;
    //VEC3D_T a4 = a*a3;
    VEC3D_T G  = 2.0*OMEGA_R*z3 + 3.0/2.0*OMEGA_M*z2 + (1.0-OMEGA_0)*z1;
    VEC3D_T H  = 6.0*OMEGA_R*z2 + 3.0*OMEGA_M*z1 + (1.0-OMEGA_0);
    VEC3D_T J  = 12.0*OMEGA_R*z1 + 3.0*OMEGA_M;
    VEC3D_T b  = -a3*G;
    VEC3D_T b2 = b*b;
    //VEC3D_T b3 = b*b2;
    VEC3D_T c  = -3.0*a2*b*G - a3*H;
    //VEC3D_T d  = -3.0*a*(2.0*b2 + a*c)*G - 6.0*a2*b*H - a3*J;

    VEC3D_T A = 1.0/a;
    VEC3D_T B = b/2.0/a3;
    VEC3D_T C = (b2/2.0/a2 - c/6.0/a)/a3;
    //VEC3D_T D = (d/24.0 - 5.0/12.0*b*c/a + 5.0/8.0*b3/a2)/a3/a2;

    VEC3D_T L_prop_s = (r_new.Norm() - Lepton->m_r4.r.Norm())/CGS_HUBRAD;
    VEC3D_T L_prop_s2 = L_prop_s*L_prop_s;
    VEC3D_T L_prop_s3 = L_prop_s*L_prop_s2;
    //VEC3D_T L_prop_s4 = L_prop_s*L_prop_s3;

    VEC3D_T delta_z_s = A*L_prop_s + B*L_prop_s2+C*L_prop_s3; //+ D*L_prop_s4;
    VEC3D_T delta_z_s2 = delta_z_s*delta_z_s;
    VEC3D_T delta_z_s3 = delta_z_s*delta_z_s2;
    //VEC3D_T delta_z_s4 = delta_z_s*delta_z_s3;

    VEC3D_T alpha = a/z1;
    VEC3D_T beta = -a/z2 + b/z1;
    VEC3D_T gamma = 2.0*a/z3 - 2.0*b/z2 + c/z1;
    //VEC3D_T delta = -6.0*a/z4 + 6.0*b/z3 - 3.0*c/z2 + d/z1;

    VEC3D_T tau_s = alpha*delta_z_s - beta/2.0*delta_z_s2 +
      gamma/6.0*delta_z_s3;// - delta/24.0*delta_z_s4;


    // Time Delay Testing:
    // std::ofstream TimeList("time_cmds.txt",std::ios::app);

//     std::string cmd_begin = "tau_strue = integrate[1/(1+z)/sqrt[OmegaR*((1+z)^4) + OmegaM*((1+z)^3) + OmegaL + (1-Omega0)*((1+z)^2)], \\ \n";
//     TimeList<<std::setprecision(32)<<cmd_begin<<"{z, "
// 	    <<(z-delta_z_s)<<", "<<z<< "}]"<<std::endl;
//     TimeList<<"tau_s = "<<tau_s<<std::endl;
//     TimeList<<"(tau_s - tau_strue)/tau_strue"<<std::endl;

//     std::cout<<(z - delta_z_s)<<" "<<z<<std::endl;

//     TimeList.close();

//     char getline;
//     std::cin>>getline;


    //------------------------ END z_s -----------------------------//



    //------------------------ BEGIN z Computation --------------------------//
    VEC3D_T delta_z = Delta_z(PL, Lepton);
    VEC3D_T delta_z2 = delta_z*delta_z;
    VEC3D_T delta_z3 = delta_z*delta_z2;
    //VEC3D_T delta_z4 = delta_z*delta_z3;

    // Compute Phi Derivatives:
    z = Lepton->m_z;
    z1 = (1.+z);
    z2 = z1*z1;
    z3 = z2*z1;
    z4 = z3*z1;
    Q  = sqrt(OMEGA_R*z4+OMEGA_M*z3+OMEGA_L+(1.0-OMEGA_0)*z2);
    VEC3D_T F = 1.0/Q;
    VEC3D_T F2 = F*F;
    VEC3D_T F3 = F*F2;
    G = 2.0*OMEGA_R*z3 + 3.0/2.0*OMEGA_M*z2 + (1.0-OMEGA_0)*z1;
    H = 6.0*OMEGA_R*z2 + 3.0*OMEGA_M*z1 + (1.0-OMEGA_0);
    J  = 12.0*OMEGA_R*z1 + 3.0*OMEGA_M;
    VEC3D_T F_pr1= -F3*G;
    VEC3D_T F_pr2= 3.0/2.0*G*F3*F2 - F3*H;
    //VEC3D_T F_pr3= -3.0*F*(2.0*F_pr1*F_pr1 + F*F_pr2)*G-6.0*F2*F_pr1*H-F3*J;
    //VEC3D_T F_pr3 = F3*(3.0/2.0*(H+J)*F2 - 15.0/4.0*G*F3*F - J );

    VEC3D_T P = (1.0 - rho/z2/2.0);
    VEC3D_T K = P/z1;
    VEC3D_T K_pr = (3.0/2.0*rho/z2 - 1.0)/z2;
    VEC3D_T K_pr2 = (2.0 - 6.0*rho/z2)/z3;
    //VEC3D_T K_pr3 = 6.0*(5.0*rho/z2 - 1.0)/z4;

    VEC3D_T Phi = K*F;
    VEC3D_T Phi_pr1 = K_pr*F + K*F_pr1;
    VEC3D_T Phi_pr2 = K_pr2*F + 2.0*K_pr*F_pr1 + K*F_pr2;
    //VEC3D_T Phi_pr3 = K_pr3*F + 3.0*K_pr2*F_pr1 + 3.0*K_pr*F_pr2 + K*F_pr3;

    VEC3D_T tau = Phi*delta_z - Phi_pr1/2.0*delta_z2 +
      Phi_pr2/6.0*delta_z3;// - Phi_pr3/24.0*delta_z4;

    // Time Delay Testing:
//     std::ofstream TimeList("time_cmds.txt",std::ios::app);

//     TimeList<<std::setprecision(32)<<"rho = "<<rho<<std::endl;

//     std::string cmd_begin = "tau_true = integrate[1/sqrt[OmegaR*((1+z)^4) + OmegaM*((1+z)^3) + OmegaL + (1-Omega0)*((1+z)^2)]/sqrt[rho+(1+z)^2], \\ \n";
//     //std::string cmd_begin = "tau_true = integrate[1/sqrt[OmegaR*((1+z)^4) + OmegaM*((1+z)^3) + OmegaL + (1-Omega0)*((1+z)^2)]/(1+z)/(1-rho/2/(1+z)^2), \\ \n";
//     TimeList<<std::setprecision(32)<<cmd_begin<<"{z, "
// 	    <<(z-delta_z)<<", "<<z<< "}]"<<std::endl;
//     TimeList<<"tau = "<<tau<<std::endl;
//     TimeList<<"(tau - tau_true)/tau_true"<<std::endl;

//     std::cout<<(z - delta_z)<<" "<<z<<" "<<delta_z<<std::endl;

//     TimeList.close();

//     char getline;
//     std::cin>>getline;

    //------------------------ END z -----------------------------//

    VEC3D_T time_delay = (tau - tau_s)/HUB_CONST;

    delta_time = time_delay;
    return delta_z_s;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  void MagneticGrid::Propagation(RelParticle& Photon,RelParticle*& Lepton,
		     VEC3D_T& PL, Vec3D& r_new,VEC3D_T& time_delay,Vec3D& n_eo,
		     VEC3D_T& delta_z,VEC3D_T& delta_zs,Vec3D& e_b)
    /*
       Updates all Lepton/Photon parameters, with these input parameters.

     \param
         Photon     - The photon to be updated
	 Lepton     - The lepton to be updated
	 r_new      - new position of Lepton/Photon
	 time_delay - time_delay of Lepton/Photon
	 n_eo       - Lepton direction before IC scat with photon.
	 delta_z    - change in redshift of the lepton over the
	              distance of the propagation length.
	 delta_zs   - change in z_s over the propagation length
	 e_b        - direction of magnetic field in current cell

    \returns void
    */
  {

    VEC3D_T ze   = Lepton->m_z;
    VEC3D_T m_eV = sqrt(Lepton->m_p4*Lepton->m_p4);
    //std::cerr<<"m_eV (Propagation): "<<m_eV<<std::endl;
    VEC3D_T p_o  = Lepton->m_p4.r.Norm()/(1.0+ze);
    //VEC3D_T rho  = m_eV*m_eV/p_o/p_o;              // [unitless]

    VEC3D_T Kappa_o = CGS_C*1.E-8;
    Kappa_o *= Lepton->m_q*m_bmag/p_o;

    VEC3D_T KappaL = Kappa_o*PL;


    //Update Lepton parameters:
    Lepton->m_r4.r = r_new;
    Lepton->m_z -= delta_z;
    Lepton->m_z_s -= delta_zs;
    Lepton->m_r4.r0 += time_delay;


    //-----------------------------------------------------------------
    // Photon rotation
    //-----------------------------------------------------------------
    Vec3D n_eo_ = (n_eo-(n_eo*e_b)*e_b)*cos(KappaL)+(n_eo^e_b)*sin(KappaL)
      +(n_eo*e_b)*e_b;

    // Unless rotation angle is very small, do rotation:
    Vec3D n_cp = n_eo_^n_eo;

    // If theta is close to Pi/2, this will avoid the domain error if arg
    //   to asin() ~ 1.0 + epsilon.
    VEC3D_T theta = "0.0";
    if (n_cp.Norm() < 0.99)
      theta = asin(n_cp.Norm());
    else {
      VEC3D_T n_dp = (n_eo_*n_eo)/n_eo_.Norm()/n_eo.Norm();
      theta = acos(n_dp);
    }

    //   VEC3D_T theta = asin(n_cp.Norm());
    if ( fabs(theta) > 1.0E-15 )
      {
	Vec3D axis = theta*n_cp/n_cp.Norm();
	Vec3D n_photon(Photon.m_p4.r/Photon.m_p4.r.Norm());
	n_photon.Rotate(axis);
	Photon.m_p4.r0 *= (1 + ze-delta_z)/(1 + ze);
	Photon.m_p4.r = Photon.m_p4.r0*n_photon;
      }
    else // No rotation, just modify 4-mom due to expansion
      {
	Photon.m_p4.r0 *= (1 + ze-delta_z)/(1 + ze);
	Photon.m_p4.r *= Photon.m_p4.r0/Photon.m_p4.r0;
      }

    Lepton->m_p4.r = p_o*(1.0 + Lepton->m_z)*n_eo_;
    Lepton->m_p4.r0 = sqrt(Lepton->m_p4.r*Lepton->m_p4.r + m_eV*m_eV);
    //std::cerr<<"Lepton 4 mom: "<<Lepton->m_p4<<std::endl<<std::endl;

    Photon.m_r4  = Lepton->m_r4;
    Photon.m_z   = Lepton->m_z;
    Photon.m_z_s = Lepton->m_z_s;


  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  ICoord MagneticGrid::CheckCurrentCell(Vec3D r)
    /*
      Calculates ICoord based on the current position in the grid. This is
      converts a VEC3D_T into a double, so it can then be converted to an
      integer to find the cell at a given particle position.

      \param  r   - position of the particle [cm]

      \return pos - ICoord, collection of 3 integers which
                    denotes the grid cell.

    */
  {

    double rx = Double(r.x);
    double ry = Double(r.y);
    double rz = Double(r.z);

    int n_x = (int) round(rx/m_cellsize);
    int n_y = (int) round(ry/m_cellsize);
    int n_z = (int) round(rz/m_cellsize);

    return ICoord(n_x,n_y,n_z);

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // CheckMagneticField
  // \param: ICoord, an index for the map
  inline Vec3D MagneticGrid::CheckMagneticField(const ICoord& c)
  {
    if(m_MagneticField.find(c) == m_MagneticField.end())
      {
        // cell has not been used before so create the field
        //B bnew = generateRandomField();
        Vec3D bnew;
        bnew = bnew.UniformSphereDirection(m_rng);
        m_MagneticField[c] = bnew;
      }

    Vec3D b = m_MagneticField[c];

    return b;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  inline Vec3D MagneticGrid::CheckMagneticField_Lock(const ICoord& c)
    /*
      Algorithm as follows:
        1) Check if ICoord defined
	2) If not, open file for reading, if ICoord in file, copy whole file
	   into MagField map, return approp B value.
	3) If ICoord not in file, open file for reading, copy whole file into
	   MagField Map, close for reading, lock immediately, open for writing,
	   define B value at this ICoord, write map to file, lock/close.
    */

  {

    if(m_MagneticField.find(c) == m_MagneticField.end())
      {

	//char* filename = sfilename.c_str();

	ICoord c_ = c;

	// Open for reading, after obtaining the lock...
	std::ifstream infile;
	int* fd = new int;
	std::string rw_type = "read";
	int error_msg = LockAttempt(m_sfilename.c_str(), fd, rw_type);
	infile.open(m_sfilename.c_str());

	if (!infile) {
	  std::cerr << "Unable to open: "<< m_sfilename <<std::endl;
	  std::cerr << "ERROR CheckMagneticField_Lock 1." << std::endl;
	  std::cerr << "error_msg = " << error_msg << std::endl;
	  exit(EXIT_FAILURE);   // call system to stop
	}

	int x = 0;
	int y = 0;
	int z = 0;
	double Bx = 0.0;
	double By = 0.0;
	double Bz = 0.0;
	ICoord coord(x,y,z);
	Vec3D Bin(Bx,By,Bz);
	while(infile) {
	  infile>>x>>y>>z;
	  infile>>Bx>>By>>Bz;

	  coord.ix = x;
	  coord.iy = y;
	  coord.iz = z;

	  // If the grid cell has been defined already in the file
	  if (c_ == coord) {
	    
	    ReadFieldMap(infile);
	    Vec3D b = m_MagneticField[c];

	    std::cerr<<"ICoord found in MF file."<<std::endl;

	    infile.close();
	    close(*fd);
	    delete fd;

	    return b;
	  }
	}
	
	// Grid cell has not been defined in file.
	// Need: close file for reading, open for writing, AND lock! Then:
	//   1) read entire map from file into MagneticField[], 
	//   2) define random field vector
	//   3) write map into text file	

	infile.clear(); // Must clear ifstream before resetting file ptr!!
	infile.seekg (0, std::ios::beg);
	infile.open(m_sfilename.c_str());

 	if (error_msg == -2){
 	  std::cerr<<"ERROR: Locking error 1."<<std::endl;
 	  exit(EXIT_FAILURE);
 	}

	ReadFieldMap(infile);

	infile.close();
	close(*fd);

	std::ofstream outfile;
 	rw_type = "write";
 	error_msg = LockAttempt(m_sfilename.c_str(), fd, rw_type);

 	if (error_msg < 0){
 	  std::cerr<<"ERROR: Locking error 2."<<std::endl;
 	  std::cerr<<"error_msg = "<<error_msg<<std::endl;
 	  exit(EXIT_FAILURE);
 	}
	
	outfile.open(m_sfilename.c_str());
	
 	if (outfile.fail()) {
 	  std::cerr << "Unable to open: "<< m_sfilename <<std::endl;
 	  std::cerr << "ERROR CheckMagneticField_Lock 2." << std::endl;
 	  std::cerr << "error_msg = "<<error_msg<<std::endl;
 	  exit(EXIT_FAILURE);
 	}

	Vec3D bnew;
	bnew = bnew.UniformSphereDirection(m_rng);
	m_MagneticField[c] = bnew;
	PrintFieldMapToFile(outfile);

	outfile.close();
	close(*fd);
	delete fd;

	return bnew;

      }
    else {
      Vec3D b = m_MagneticField[c];
      return b;
    }

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //  void MagneticGrid::ReadFieldMap(Field& MagneticField, 
  void MagneticGrid::ReadFieldMap( std::ifstream& infile)
    /*
      
    Reads data from filename textfile into m_MagneticField map.

    */
  {

    infile.clear();
    infile.seekg (0, std::ios::beg);

    if (!infile) {
      std::cerr << "Unable to open infile." <<std::endl;
      std::cerr << "ERROR ReadFieldMap." << std::endl;
      exit(EXIT_FAILURE);   // call system to stop
    }

    int x = 0;
    int y = 0;
    int z = 0;
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;
    ICoord coord(x,y,z);
    Vec3D Bin(Bx,By,Bz);
    while(infile){

      infile>>x>>y>>z;
      infile>>Bx>>By>>Bz;
    
      coord.ix = x;
      coord.iy = y;
      coord.iz = z;

      Bin.x = Bx;
      Bin.y = By;
      Bin.z = Bz;
      m_MagneticField[coord] = Bin;

    }

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  void MagneticGrid::PrintFieldMapToFile(std::ofstream& outfile)
    /*
      Writes out m_MagneticField Map data into textfile filename.
    */
  {

    if (!outfile) {
      std::cerr<<"Unable to open outfile."<<std::endl;
      std::cerr << "ERROR PrintFieldMapToFile." << std::endl;
      exit(EXIT_FAILURE);
    }
  
    for(Field::const_iterator it = m_MagneticField.begin();
	it != m_MagneticField.end(); ++it)
      {
	ICoord coord = it->first;
	Vec3D bvalue = it->second;
	outfile << coord <<" ";
	outfile << bvalue.x<<" "<<bvalue.y<<" "<<bvalue.z << std::endl;
      }

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  int MagneticGrid::LockAttempt(const char* filename, int* fd,
				std::string rwtype)
    /*
      file locking (so multiple users don't interfere with each other)
      
      \return: a file descriptor if successful, else negative number for
               error message.
      
      \note: locks are released when the files are closed in UnLock()
             function when file descriptor is closed.
    */
  {

    const int WAIT_TIME = 10;        // 10 musec
    struct flock lck;

    if (rwtype == "read") {  // read_write = TRUE for read lock
      *fd = open(filename, O_RDONLY);
      lck.l_type = F_RDLCK;
      std::cerr<<"Read lock attempted."<<std::endl;
    } else if (rwtype == "write") {
      *fd = open(filename, O_WRONLY);
      lck.l_type = F_WRLCK;
      std::cerr<<"Write lock attempted."<<std::endl;
    } else {
      std::cerr<<"Invalid Lock Type: "<<rwtype<<std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (*fd < 0) return -1; // file not found, or could not be opened
    
    // set up the record locking structure, the address of which
    // is passed to the fcntl system call.
    lck.l_whence = SEEK_SET;
    lck.l_start = 0;                            // from beginning of file
    lck.l_len = 0;                              // until end of the file
    
    // Attempt locking

    ///////////////////////////////////////////////////////////////
    // NOTE: errno codes are listed in /usr/include/errno.h
    //   Check them to see what the problem was if needed! (or online)
    ///////////////////////////////////////////////////////////////

    int TRYNO = 0;
    std::cerr<<"Lock Attempt..."<<std::endl;
    while ( TRYNO >= 0 ) { // only a return can break out of the while loop

      if ( fcntl(*fd, F_SETLK, &lck) == 0 ){
	std::cerr<<"TRYNO: "<<TRYNO<<std::endl;
	return *fd;
      }
      else if (errno == EAGAIN || errno == EACCES) { // file in use
	//std::cerr<<"File in use..."<<std::endl;
	usleep(WAIT_TIME);
      }
      else if (errno == EIO) {
	std::cerr<<"IO Locking ERROR? errno: "<<errno<<std::endl;
	usleep(WAIT_TIME);
      }
      else {
	std::cerr<<"Locking ERROR? errno: "<<std::endl;
	std::cerr<<"TRYNO: "<<TRYNO<<std::endl;
	usleep(WAIT_TIME);
	//std::cerr<<"errno: "<<errno<<std::endl;
	//return -2;
      }
      TRYNO = TRYNO+1;

    }

    return -3; // nothing happened!

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


}
