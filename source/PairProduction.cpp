/*!
-------------------------------------------------------------------------------
    \file   PairProduction.cpp

     Pair Production class implementation file for FRW Universe.

    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
         arlen@astro.ucla.edu                  \n

    \date      July 4, 2007


    \revision: 02/15/2008 - Reorganized constants
               07/09/2008 - Changed PropagationLengthEBL() to VEC3D_T
                      inputs instead of double
         08/13/2008 - Changed functions to input RelParticle instead
                      of their components.
         08/19/2008 - Changed (almost) everything to VEC3D_T, especially
                      affects PropagationLengthEBL().
          Only functions we keep with double precision are
          in DIRBR.cpp.
         10/17/2008 - Rewrite the PropagationLength() function so that
                      it "stops" integrating when Photon has reached
          R_p (radius in comoving coords).
         6/7/2009   - Changed constructor to initialize DIRBR class
                      only once per run.
         4/14/2010  - Changed the loop condition on the photon integrate
                      to next z pair production point. Needed to add the
          map to keep track of the exact distance as well.

         10/7/2010  - Removed all dependency on RelParticle class.
                      Should be general enough to work with all kinds of
          inputs (not relying on RelParticle).

    \note
            In this version, for testing the geometry of the photon, I needed
      to modify the while loop for exiting when L <= PropLength, instead
      of R_proper.
-------------------------------------------------------------------------------
*/

#include "PairProduction.hpp"
#include <iomanip>
#include <map>

using PhysConst::OMEGA_0;
using PhysConst::OMEGA_L;
using PhysConst::OMEGA_M;
using PhysConst::OMEGA_R;

namespace IGCascade {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _rng: rng
PairProduction::PairProduction(RandomNumbers *_rng, VEC3D_T ze) {

  m_rng = _rng;

  // We need type double versions of these, because they occur in
  // the integration loop of computing tau for speed.
  m_OmegaR_d = Double(OMEGA_R);
  m_OmegaM_d = Double(OMEGA_M);
  m_OmegaL_d = Double(OMEGA_L);
  m_Omega0_d = Double(OMEGA_0);

  DefineNumericConst();
  DefineIntegrationParameters(ze);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void PairProduction::DefineNumericConst(void)
/*!
  Defines numeric constants used in PairProduction class. Only called from
  Constructor.
*/
{

  // Public member:
  m_DE = "1.0E-25";

  //~~~~~~~~VEC3D_T Values~~~~~~~~~
  D0 = "0.0";
  D1 = "1.0";
  D2 = "2.0";
  D3 = "3.0";
  D4 = "4.0";
  D5 = "5.0";
  D6 = "6.0";
  D7 = "7.0";
  D8 = "8.0";

  m_DE_int = "1.0E-12";                    // Rel comp precision of roots
  m_PI2d6 = VEC3D_PI * VEC3D_PI / D2 / D3; // pi*pi/6 for PolyLog func
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void PairProduction::DefineIntegrationParameters(VEC3D_T ze) {

  /////////////////////////////////////////////////
  // If we are very close to the source, we'd like more steps to integrate
  // over to find the interaction distance.
  if (ze < 5.0e-2) m_dz_dd = ze / 50.0; // 50 integration steps.
  // if (ze < 2.0e-2) m_dz_dd = "1.0e-4";
  else {                // propagation length dz
    m_dz_dd = "1.0e-3"; // ~3 Mpc dist
    // m_dz_dd = "3.0e-4";     // ~1 Mpc dist
    // m_dz_dd = "1.0e-4";
  }
  VEC3D_T scale = "10.0";
  m_dz_time_delay = m_dz_dd / scale; // time_delay dz
  // m_dz_time_delay = VEC3D_T("1.223961751575643313e-05");

  //~~~~~~~~~double values~~~~~~~~~~~~
  // Integration parameters for innermost PropagationLengthEBL integral:
  double qmin = 5.0e-4;
  m_dlog_q = 2.0e-4;
  m_dz = Double(m_dz_dd);
  m_trans_prob_min = 0.1; // 10 % survival probability.

  m_imax = -(int)(log(qmin) / m_dlog_q);
  // NOTE: m_imax = 38004

  double f_i = 0.0;
  double yi = 0.0;
  // Set up logarithmic grid for the F(q) integration.
  for (int i = 0; i < m_imax; i++) {
    yi = exp(-(double)(i)*m_dlog_q);
    m_q_vec.push_back(yi); // q[0]=1.0;

    f_i += ((1.0 + yi - yi * yi / 2.0) *
                log((1.0 + sqrt(1.0 - yi)) / (1.0 - sqrt(1.0 - yi))) -
            (1.0 + yi) * sqrt(1.0 - yi)) /
           yi;
    m_F_vec.push_back(2.0 * m_q_vec[i] * m_q_vec[i] * f_i * m_dlog_q);
  }
  //////////////////////////////////////////////////////////////
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool PairProduction::CheckPairProductionEBL(
    DIRBRBase *ebl_model, const double gam_ph_egy, const VEC3D_T z_o,
    const double z_min, VEC3D_T &z_int, double &TotalLambdaInt
)
/*!
  Checks if pair production occurs through ebl_model with a gamma photon of
  gam_ph_egy and starting redshift of z_o, by throwing a random number.
  Calls IntegrateToTau() to do the heavy lifting.

  MODIFIES:
    1) z_int - the interaction redshift
    2) TotalLambdaInt - the total integral over ebl photon wavelength, lambda
       which is used later in determining with what energy of ebl photon
 the Gamma Photon pair produced with.

  \returns - bool true or false depending on whether pair production took
             place.
 */
{

  double chi = m_rng->Uniform();
  double lnchi = log(chi);
  const double tauFinal = -lnchi;

  return IntegrateToTau(
      ebl_model, gam_ph_egy, z_o, z_min, tauFinal, z_int, TotalLambdaInt
  );
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// void PairProduction::
// GetOpticalDepth(DIRBRBase* ebl_model,const double gam_ph_egy,
//		  const VEC3D_T z_o, const VEC3D_T z_final)
//{

//}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool PairProduction::IntegrateToTau(
    DIRBRBase *ebl_model, const double gam_ph_egy, const VEC3D_T z_o,
    const double z_min, const double tauFinal, VEC3D_T &z_int,
    double &TotalLambdaInt
)
/*
  Function that does the heavy lifting of the tau (Optical Depth)
  integration.

  NOTE:
      1) Uses main while loop in double precision (rather than double-double
      to speed things up and accuracy here is not needed) to determine approx.
  prop_steps (<=> change in redshift) to where photon interacted. The dz
  step here is m_dz (of the class) and its precision depends on the
  redshift.

2) The final step interpolates the redshift step from the
crude delta_z in the integration to a delta_z between these
bins. (ADDED 4/16/2013)

*/
{

  /////////////////////////////////////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~Propagation Length Integral~~~~~~~~~~~~~~~~~~~~//
  // NOTE: Cannot use the convert to double function, Double()       //
  //   within these loops because it slows the computation way down  //
  //   So make new double variables at start of this subroutine.     //
  /////////////////////////////////////////////////////////////////////
  // NOTE: [tau_const] = [m^2 s sr eV^-1 ], Integration Prefactor
  double hc_d = Double(PhysConst::eVCM_HC * 1.0E4); // [eV cm] -> [eV mum]
  double me_eV_d = Double(PhysConst::eV_MELEC);
  double me_sq = me_eV_d * me_eV_d;
  double nJtoeV_d = Double(PhysConst::nJ_TO_eV);
  double tau_const = Double(
      3.0 * PhysConst::CGS_THOM_CS * 1.0E-4 * gam_ph_egy * 4.0 * VEC3D_PI /
      (8.0 * PhysConst::HUB_CONST * me_eV_d * me_eV_d * (1.0 + z_o))
  );
  double DimlessTotalLambdaInt = 0.0;
  TotalLambdaInt = 0.0;
  double Q = 0.0;
  double tau = 0.0;
  double z_o_d = Double(z_o);
  double z = z_o_d;

  // Find Interaction point (low precision)
  double tau_prior = 0.0;
  double z_prior = 0.0;
  // VEC3D_T prop_steps = 0.0;
  while (tau < tauFinal) { // Also exits when z < z_min (below)
    double z1 = (1. + z);
    double z2 = z1 * z1;
    double z3 = z1 * z2;
    double z4 = z1 * z3;
    Q = sqrt(
        m_OmegaR_d * z4 + m_OmegaM_d * z3 + m_OmegaL_d + (1.0 - m_Omega0_d) * z2
    );

    double dirbr_temp = 1.0;
    DimlessTotalLambdaInt = 0.0;

    // F(q) and DIRBR integral
    double lambda_const = hc_d * gam_ph_egy * z2 / (me_sq * (1.0 + z_o_d));
    for (int i = (m_imax - 1); i >= 0; i--) {
      double lambda = 0.0;
      lambda = lambda_const * m_q_vec[i];
      dirbr_temp = ebl_model->GetDIRBR(lambda, z);
      if (lambda > 100.0) dirbr_temp += ebl_model->GetCMBR(lambda);
      DimlessTotalLambdaInt += m_q_vec[i] * m_F_vec[i] * dirbr_temp;
    }
    TotalLambdaInt = DimlessTotalLambdaInt * m_dlog_q * nJtoeV_d;
    // tau+= TotalLambdaInt*z4*tau_const*m_dz/Q;
    tau_prior = tau;
    tau = tau_prior + TotalLambdaInt * z4 * tau_const * m_dz / Q;
    z_prior = z;

    z -= m_dz;
    // prop_steps += 1.0;

    if (z < z_min) {
      z_int = z_min;
      return false;
    }
  }
  /////////////////////////////////////////////////////////////////////

  // Perform interpolation:
  z_int = (VEC3D_T)(z_prior +
                    (z - z_prior) * (tauFinal - tau_prior) / (tau - tau_prior));

  // interacted redshift in double-double precision.
  // z_int = (z_o - m_dz*prop_steps);
  std::cout << "\n  z_int: " << z_int << std::endl << std::endl;
  return true;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool PairProduction::UpdateGammaPhoton(
    Vec4D &gam_ph_p4, Vec4D &gam_ph_r4, VEC3D_T &gam_ph_z, VEC3D_T &gam_ph_z_s,
    VEC3D_T &delta_z_step
)
/*!  Updates Gamma's .r4, .z, .z_s, and .p4

  \param
      gam_ph_p4  - gamma ray 4 momentum
gam_ph_r4  - gamma ray 4 position
gam_ph_z   - gamma ray redshift
gam_ph_z_s - gamma ray z_s, ( ~ time delay; see Cosmology writeup!)
delta_z_step - redshift difference that photon took on this step.

   \return True - if pair production occurs before z=0 False otherwise.
*/
{

  bool pair_prod_occurs = true;

  // VEC3D_T z_o = gam_ph_z;
  VEC3D_T delta_z = m_dz_time_delay;

  Vec3D R_o = gam_ph_r4.r;                       // Photon initial radial coord
  Vec3D e_ph = gam_ph_p4.r / gam_ph_p4.r.Norm(); // Photon direction
  VEC3D_T z_final = gam_ph_z - delta_z_step;

  //////////////////////////////////////////////////////////////
  //------------------Time Delay Computation------------------//
  //////////////////////////////////////////////////////////////
  // VEC3D_T Delta_time = "0.0";
  VEC3D_T L_prop = "0.0"; // Dist traveled during propagation
  VEC3D_T z_dd = gam_ph_z;
  unsigned nsteps = 0;
  // while ( ((z_dd-z_final) > m_DE) && (gam_ph_z_s > m_DE) ) {
  while ((z_dd > z_final) && (gam_ph_z_s > m_DE)) {

    if ((z_dd - z_final) < delta_z) delta_z = (z_dd - z_final);

    // std::cout<<"  z_dd: "<<z_dd<<" z_final: "<<z_final<<" delta_z:
    // "<<delta_z<<std::endl;
    VEC3D_T Delta_time = "0.0";
    VEC3D_T deltaz_s = "0.0";

    NonRadialPhotonPropagation(
        gam_ph_p4, gam_ph_r4, gam_ph_z, gam_ph_z_s, delta_z, L_prop, deltaz_s,
        Delta_time
    );

    // If z_s = 0 surface crossed:
    if ((gam_ph_z_s - deltaz_s) < -m_DE) {
      // std::cout<<"z_s = 0.0 surface crossed."<<std::endl;
      VEC3D_T delta_z_L = "0.0";
      VEC3D_T delta_z_R = delta_z;

      while (fabs(gam_ph_z_s - deltaz_s) > m_DE) {

        delta_z = (delta_z_L + delta_z_R) / 2.0;

        NonRadialPhotonPropagation(
            gam_ph_p4, gam_ph_r4, gam_ph_z, gam_ph_z_s, delta_z, L_prop,
            deltaz_s, Delta_time
        );

        if ((gam_ph_z_s - deltaz_s) > 0.0)
          delta_z_L = delta_z;
        else
          delta_z_R = delta_z;
      }

      pair_prod_occurs = false;
    }

    gam_ph_z_s -= deltaz_s;
    gam_ph_r4.r0 += Delta_time;
    gam_ph_r4.r += L_prop * PhysConst::CGS_HUBRAD * e_ph;

    gam_ph_p4.r0 *= (1.0 + z_dd - delta_z) / (1.0 + z_dd);
    gam_ph_p4.r *= gam_ph_p4.r0 / gam_ph_p4.r.Norm();

    z_dd -= delta_z;
    gam_ph_z = z_dd;

    nsteps++;
  }

  return pair_prod_occurs;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VEC3D_T PairProduction::GetEBLPhotonEgy(
    DIRBRBase *ebl_model, const double TotalLambdaInt, const double z_int,
    const double z_o, const double gam_ph_egy
)
/*!
  Main purpose is to define the energy of the EBL Photon which interacts
  with the energetic Gamma Photon, by performing innermost integral of
  tau integration.
  Results in updated background photon's energy.

  \param
       z_int          - redshift at which gamma photon interacts.
 z_o            - redshift at which gamma photon started.
 gam_ph_egy     - gamma photon energy
 bg_ph_egy      - background photon energy
 TotalLambdaInt - total lambda integral, as calculated in
                  CheckPairProductionEBL() fn.

  \returns - ebl photon egy which interacted to pair produce
*/
{

  VEC3D_T bg_ph_egy = "0.0";

  ///////////////////////////////////////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~Energy of EBL Photon Integral~~~~~~~~~~~~~~~~~~//
  // NOTE: This determines the EBL Photon energy at z=0, not where it  //
  //   interacts with gamma photon. Thus, scale accordingly at the end //
  ///////////////////////////////////////////////////////////////////////
  double z2 = (1. + z_int) * (1. + z_int);
  double hc_d = Double(PhysConst::eVCM_HC * 1.0E4); // [eV cm] -> [eV mum]
  double me_eV_d = Double(PhysConst::eV_MELEC);
  double me_sq = me_eV_d * me_eV_d;
  double nJtoeV_d = Double(PhysConst::nJ_TO_eV);
  double LambdaInt = 0.0;
  double EBL_lambda = 0.0;
  double EBL_lambda_next = 0.0;
  double EPS_dirbr = 1.0e-10;

  double chi_rand = m_rng->Uniform();

  int i = 0;
  LambdaInt = 0.0;
  EBL_lambda = m_q_vec[i] * hc_d * gam_ph_egy * z2 / (me_sq * (1.0 + z_o));

  double dirbr_temp = ebl_model->GetDIRBR(EBL_lambda, z_int);
  if (EBL_lambda > 100.0) dirbr_temp += ebl_model->GetCMBR(EBL_lambda);

  // double dirbr = nJtoeV_d*TimeEvolDIRBR(EBL_lambda,z,dirbr_temp);
  //  Note: units are in [ eV m^-2 sr^-1 s^-1 ]
  double dirbr = nJtoeV_d * dirbr_temp;

  std::cout << std::endl << "Integrating to find EBL Photon energy...\n\n";

  while (LambdaInt < TotalLambdaInt * chi_rand) {

    EBL_lambda_next =
        m_q_vec[i + 1] * hc_d * gam_ph_egy * z2 / (me_sq * (1.0 + z_o));

    double dirbr_temp_next = ebl_model->GetDIRBR(EBL_lambda_next, z_int);
    if (EBL_lambda_next > 100.0) {
      dirbr_temp_next += ebl_model->GetCMBR(EBL_lambda_next);
    }

    // double dirbr_next = nJtoeV_d*TimeEvolDIRBR(EBL_lambda_next,z,
    //       dirbr_temp_next);       // [ eV m^-2 sr^-1 s^-1 ]
    double dirbr_next = nJtoeV_d * dirbr_temp_next; // [ eV m^-2 sr^-1 s^-1 ]

    // Interpolation of m_q_vec[i]
    if (dirbr_next <= EPS_dirbr) {
      LambdaInt = 2.0 * TotalLambdaInt;
      std::cout << std::endl
                << "In between m_q_vec[i] bins"
                << "\t m_q_vec[i] before adjustment = " << m_q_vec[i]
                << "\tchi_rand = " << chi_rand << std::endl;
      m_q_vec[i] = (m_q_vec[i + 1] + m_q_vec[i]) / 2.0;
      std::cout << "m_q_vec[i] after adjustment = " << m_q_vec[i] << std::endl
                << std::endl;
    } else {
      LambdaInt +=
          ((m_q_vec[i] * m_F_vec[i] * dirbr +
            m_q_vec[i + 1] * m_F_vec[i + 1] * dirbr_next) *
           m_dlog_q / 2.0);
      i++;

      EBL_lambda = EBL_lambda_next;
      dirbr_temp = dirbr_temp_next;
      dirbr = dirbr_next;
    }
  }
  //////////////////////////////////////////////////////////////

  bg_ph_egy = (PhysConst::eV_MELEC * PhysConst::eV_MELEC * (1.0 + z_o)) /
              gam_ph_egy / m_q_vec[i - 1] / z2;
  bg_ph_egy *= (1.0 + z_int);

  return bg_ph_egy;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void PairProduction::UpdateEBLPhoton(
    Vec4D &gam_ph_p4, Vec4D &gam_ph_r4, VEC3D_T &gam_ph_z, VEC3D_T &gam_ph_z_s,
    Vec4D &bg_ph_p4, Vec4D &bg_ph_r4, VEC3D_T &bg_ph_z
)
/*! Routine defines the EBL 3 momentum, by calling the ImpactAngle()
  function and rotating into the frame where direction of gamma
  photon is along z_axis.
  Also updates all of EBLPhoton dynamical parameters.

  \param
      gam_ph_p4  - gamma ray 4 momentum
gam_ph_r4  - gamma ray 4 position
gam_ph_z   - gamma ray redshift
gam_ph_z_s - gamma ray z_s, ( ~ time delay; see Cosmology writeup!)
bg_ph_p4   - background photon 4 momentum; calculated here
bg_ph_r4   - background photon 4 position; calculated here
bg_ph_z    - background photon redshift; calculated here
*/
{

  bg_ph_z = gam_ph_z;
  bg_ph_r4.r0 = gam_ph_r4.r0;
  bg_ph_r4.r = gam_ph_r4.r;

  ////////////////////////////////////////////////////////////////
  /// Now for the hard part of updating EBL Photon 3-momentum: ///
  ////////////////////////////////////////////////////////////////
  VEC3D_T q =
      PhysConst::eV_MELEC * PhysConst::eV_MELEC / bg_ph_p4.r0 / gam_ph_p4.r0;
  if (q > D1) {
    std::cerr << "ERROR: Invalid value of q = " << q << std::endl;
    exit(EXIT_FAILURE);
  }
  // Generate polar, azimuthal angles.
  VEC3D_T cos_theta = ImpactAngle(q);
  VEC3D_T sin_theta = sqrt(D1 - cos_theta * cos_theta);
  VEC3D_T phi = (D2 * VEC3D_PI) * ((VEC3D_T)m_rng->Uniform());
  VEC3D_T sin_phi = sin(phi);
  VEC3D_T cos_phi = cos(phi);

  Vec3D n_gp = gam_ph_p4.r / gam_ph_p4.r.Norm();
  Vec3D n_ebl(D0, D0, D0);
  if (n_gp.z > D0) {
    n_ebl.x = sin_theta * cos_phi;
    n_ebl.y = sin_theta * sin_phi;
    n_ebl.z = cos_theta;
  } else {
    n_ebl.x = sin_theta * cos_phi;
    n_ebl.y = -sin_theta * sin_phi;
    n_ebl.z = -cos_theta;
  }

  ////////////////////////////////////////////////////////////
  // Unless gamma propagates parallel to z-axis (almost never)
  // rotate ebl from frame where gamma propagates along z.
  Vec3D e_z(D0, D0, D1);
  VEC3D_T CheckParallel = e_z * n_gp;
  VEC3D_T eps = 1.E-25;
  if (fabs(CheckParallel) < (D1 - eps)) {
    VEC3D_T e_xy = sqrt(n_gp.x * n_gp.x + n_gp.y * n_gp.y);
    VEC3D_T psi = atan(e_xy / n_gp.z);
    Vec3D axis(n_gp.y, -n_gp.x, D0);
    axis *= psi / e_xy;
    n_ebl.Rotate(-axis);
  }

  bg_ph_p4.r = n_ebl * bg_ph_p4.r0 / n_ebl.Norm();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VEC3D_T PairProduction::ImpactAngle(VEC3D_T q)
/*! This routine Samples the impact angle between the incident gamma
ray and background photon in the lab reference frame.

\param
  q = m^2/(E1*E2) Where E1, E2 are the gamma and background
      photon energy in the lab r.f.    [1]

  q ~ 1 - "soft" photon-near threshold energy
  q ~ 0 - "hard" photon. (gamma energy in lab frame -> infinity)

\return  cos(theta) where theta is the polar scattering angle
       between each particle's momenta in the lab r.f.  [1]

*/
{

  // Factors in the G(q) function:
  VEC3D_T phi = sqrt(D1 - q);
  VEC3D_T etaminus = D0;

  if (q < 0.001)
    etaminus = q / D2 / (D1 + phi);
  else
    etaminus = (D1 - phi) / D2;

  VEC3D_T etaplus = D1 - etaminus;
  VEC3D_T g1 = (q / D2 - D1 + D1 / q);
  VEC3D_T g2 = phi * (D2 / q - D1);
  VEC3D_T z = etaplus / etaminus;
  VEC3D_T zinv = etaminus / etaplus;

  VEC3D_T G =
      (g1 * log(z) - g2 + PolyLog1(zinv) - PolyLog1(z) + log(zinv) * log(q / D4)
      );

  VEC3D_T DEps = (VEC3D_T("1.E-12")) * G; // Gives us 12 s.f. accuracy

  VEC3D_T chi = (VEC3D_T)(m_rng->Uniform());

  VEC3D_T LB = q;
  VEC3D_T RB = D1;
  VEC3D_T eps = DEps + DEps;
  VEC3D_T x = RB;
  VEC3D_T F = D0;

  VEC3D_T etaminus_ = D0;
  VEC3D_T etaplus_ = D0;

  while (fabs(eps) > DEps) {
    x = (LB + RB) / D2;

    if (x < 0.001)
      etaminus_ = x / D2 / (D1 + sqrt(D1 - x));
    else
      etaminus_ = (D1 - sqrt(D1 - x)) / D2;

    etaplus_ = D1 - etaminus_;

    F = (x / D2 - D1 + D1 / x) * log(etaplus_ / etaminus_) -
        sqrt(D1 - x) * (D2 / x - D1) + PolyLog1(etaminus_ / etaplus_) -
        PolyLog1(etaplus_ / etaminus_) +
        log(etaminus_ / etaplus_) * log(x / D4);

    eps = F - G * chi;

    if (eps > D0)
      LB = x;
    else
      RB = x;
  }

  VEC3D_T cos_theta = D1 - D2 * q / x;

  return cos_theta;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool PairProduction::RelativisticKinematics(
    Vec4D &gam_ph_p4, Vec4D &bg_ph_p4, Vec4D &elec_p4, Vec4D &pos_p4
)
/*! Computes relativistic kinematics of pair production of electron
  and positron by two incident photons. Scattering angle is
  sampled utilizing full qed cross-section.

  \param
      gam_ph_p4  - gamma ray 4 momentum
bg_ph_p4   - background photon 4 momentum
elec_p4    - electron 4 momentum
pos_p4     - positron 4 momentum

  \note All input parameters are assumed to be in the lab frame.

*/
{

  std::cout << "Computing relativistic kinematics of interaction..."
            << std::endl;

  VEC3D_T me_sq = (PhysConst::eV_MELEC) * (PhysConst::eV_MELEC); // [eV^2]

  // Compute boost for the CM frame:
  Vec3D beta = (gam_ph_p4.r + bg_ph_p4.r) / (gam_ph_p4.r0 + bg_ph_p4.r0);

  // Error check
  if (beta.Norm() >= D1) {
    std::cerr << std::endl << "ERROR: beta must be less than 1.0" << std::endl;
    exit(EXIT_FAILURE);
  }

  VEC3D_T betasq = beta * beta;
  VEC3D_T gamma = D1 / sqrt(D1 - betasq);
  Vec3D e1 = gam_ph_p4.GetDirection();
  Vec3D e2 = bg_ph_p4.GetDirection();

  VEC3D_T Ecm =
      (gamma * gam_ph_p4.r0 * bg_ph_p4.r0 / (gam_ph_p4.r0 + bg_ph_p4.r0)) *
      (D1 - e1 * e2);
  // VEC3D_T y = D0;

  // First check if enough energy for interaction to occur:
  if (Ecm < PhysConst::eV_MELEC) {
    std::cerr << std::endl;
    std::cerr << "ERROR: Not enough energy for pair production" << std::endl;
    std::cerr << "Photon Energy - mass electron = " << Ecm - PhysConst::eV_MELEC
              << " eV";
    std::cerr << std::endl << "Ecm = " << Ecm << std::endl;
    return false;
  } else {

    // Set up coords for collision in CM frame:
    Vec3D n1 = gam_ph_p4.r0 / Ecm *
               (e1 - (beta * e1) * beta / betasq +
                gamma * (beta * e1 / betasq - D1) * beta);

    VEC3D_T betan1 = beta * n1;
    // VEC3D_T betan1sq = betan1*betan1;
    Vec3D v = beta ^ n1;
    VEC3D_T VecCompare = v.Norm() / beta.Norm();

    Vec3D n2(D1, D0, D0);

    // Set up coords in CM frame. Check to see if beta is parallel to n1
    if (fabs(VecCompare) > 1.E-12) {
      Vec3D cp = e1 ^ beta;
      VEC3D_T denom = gam_ph_p4.r0 / Ecm * cp.Norm();
      n2 = (beta - betan1 * n1) / denom;
    } else {
      Vec3D n_2(D1, D0, D0);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Probably overkill, but this makes absolutely certain that
      // we can get an orthogonal vector.
      VEC3D_T DEps = 1.E-10;
      while (fabs(VecCompare) < DEps) {
        std::cout << "v not perpendicular to n1. Retrying..." << std::endl;
        n2 = n_2.UniformSphereDirection(m_rng);
        v = n2 ^ n1;
        VecCompare = v.Norm();
      }
      n2 = (n2 - (n2 * n1) * n1) /
           sqrt(n2 * n2 - (n2 * n1) * (n2 * n1)); // normalize
    }

    Vec3D n3 = n1 ^ n2;

    VEC3D_T s = (gam_ph_p4 + bg_ph_p4) * (gam_ph_p4 + bg_ph_p4);
    VEC3D_T x = D4 * me_sq / s;

    // Azimuthal scattering angle in CM frame
    VEC3D_T phi = (D2 * VEC3D_PI) * ((VEC3D_T)m_rng->Uniform());

    // Determine polar angle of outgoing particle in CM frame
    // std::cout<<"  Get Scattering angle: "<<std::endl;
    VEC3D_T cos_theta = Scattering(x);
    VEC3D_T sin_theta = sqrt(D1 - cos_theta * cos_theta);
    // std::cout<<"  Finished."<<std::endl;

    VEC3D_T sin_phi = sin(phi);
    VEC3D_T cos_phi = cos(phi);

    VEC3D_T p_e = sqrt(Ecm * Ecm - me_sq); // Mag of 3 mom outgoing elec/pos
    Vec3D n_e = cos_theta * n1 + sin_theta * (cos_phi * n2 + sin_phi * n3);
    Vec3D pe = p_e * n_e; // 3 momentum of electron in CM frame
    Vec3D ppos = -pe;

    // Reverse transformation to get elec & pos 4-mom into lab frame:
    elec_p4.r0 = gamma * (Ecm + beta * pe);
    elec_p4.r = pe - (pe * beta) * beta / betasq +
                gamma * (pe * beta / betasq + Ecm) * beta;
    pos_p4.r0 = gamma * (Ecm + beta * ppos);
    pos_p4.r = ppos - (ppos * beta) * beta / betasq +
               gamma * (ppos * beta / betasq + Ecm) * beta;

    /////////////////////////////////////////////////////
    //--------Testing for Momentum Conservation--------//
    //---------Test precision to 12 sig. figs.---------//
    /////////////////////////////////////////////////////
    VEC3D_T s_before = (gam_ph_p4 + bg_ph_p4) * (gam_ph_p4 + bg_ph_p4);
    VEC3D_T s_after = (elec_p4 + pos_p4) * (elec_p4 + pos_p4);
    VEC3D_T eps = 1.E-12;
    if (fabs(s_before - s_after) > eps * s_after) {
      std::cout << std::setprecision(12)
                << " WARNING! Mandelstam variable s NOT conserved! ";
      std::cout << std::endl << "s_before = " << s_before << std::endl;
      std::cout << "s_after  = " << s_after << std::endl;
      std::cout << "Energy of Gamma Ray = " << gam_ph_p4.r0 << std::endl;
    }

    std::cout << "Finished relativistic kinematics." << std::endl;
    return true;

  } // end if
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VEC3D_T PairProduction::Scattering(VEC3D_T x)
/*! Samples Scattering angle of the electron (and thus the positron)
  in the particle's r.f.

 \param x = 4m^2/s. s is Mandelstam invariant = 4*E^2 where E is the
        photon energy in the CM frame.
        x -> 1 "soft" photon-near threshold energy
        x -> 0 "hard" photon. E (gamma energy in lab frame) -> infinity

 /return y = ( 1 - x )^(1/2)cos(theta) where theta is
    polar scattering angle of the electron in the CM frame
*/
{

  // Factors in the f functions
  VEC3D_T lambda = sqrt(D1 - x);
  VEC3D_T h1 = D0;
  VEC3D_T h2 = lambda / (D1 + x - x * x / D2);
  VEC3D_T j1 = D0;
  VEC3D_T j2 = D0;
  VEC3D_T logx = log(x);
  VEC3D_T sqrtx = sqrt(x);
  VEC3D_T loglam = log(D1 + lambda);

  VEC3D_T chi = (VEC3D_T)(m_rng->Uniform());

  // Chi dependent Prefactors of the f functions.
  h1 = (sqrtx / (D1 + lambda)) * exp(D2 * chi * loglam) / exp(chi * logx);
  j2 = chi * (D1 + x);

  VEC3D_T eps = D1;
  VEC3D_T LB = -D1;
  VEC3D_T RB = D1;
  VEC3D_T y = (LB + RB) / D2;

  VEC3D_T f1 = D0;
  VEC3D_T f2 = D0;

  while (fabs(eps) > m_DE_int) {

    y = (LB + RB) / D2;

    if (y < D0)
      j1 = D1 - lambda * y;
    else
      j1 = ((D1 - y) * (D1 + y) + x * y * y) / (D1 + lambda * y);

    f1 = h1 * j1 / sqrt((D1 - y) * (D1 + y) + x * y * y);
    f2 =
        exp(h2 * (j2 - (y + D1) / D2 *
                           (x * (D1 - y + x * y) /
                                (x * y * y + (D1 - y) * (D1 + y)) +
                            D1)));

    eps = f1 - f2;

    if (eps > D0)
      LB = y;
    else
      RB = y;

  } // end while

  return y;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VEC3D_T PairProduction::PolyLog1(VEC3D_T x)
/*! Dilogarithm related function. In terms of traditional notation for
    dilogarithm, Li_2(z)=z+z^2/4+z^3/9+...+z^n/n^2, it is defined as

         PolyLog1(x)=-Li_2(-x)

  In terms of integral representation

           PolyLog1(x)=Int_1^x  ln(1+x)dx/x+pi^2/12

  The function is defined for x from 0 to infinity and PolyLog1(1)=pi^2/12

*/
{
  if (x < D0) {
    std::cout << "PolyLog1: Negative argument." << std::endl;
    exit(0);
    // VEC3D_T PL=D0;
    // return PL;
  }

  if (x > D2) { // large x approximation  (x>2)

    VEC3D_T Ln = log(x);
    VEC3D_T xPower = D1;
    VEC3D_T PL = Ln * Ln / D2 + m_PI2d6;
    VEC3D_T dPL = m_PI2d6;
    VEC3D_T i = D0;

    while (fabs(dPL) > m_DE_int * m_PI2d6) {
      i += D1;
      xPower /= (-x);
      dPL = xPower / i / i;
      PL += dPL;
    }
    // std::cout<<" 1 "<<PL<<std::endl;
    return PL;

  } else if (x < ((VEC3D_T)0.618)) { // small x approximation

    VEC3D_T xPower = -D1;
    VEC3D_T PL = D0;
    VEC3D_T dPL = m_PI2d6;
    VEC3D_T i = D0;

    while (fabs(dPL) > m_DE_int * m_PI2d6) {
      i += D1;
      xPower *= (-x);
      dPL = xPower / i / i;
      PL += dPL;
    }
    // std::cout<<" 2 "<<PL<<std::endl;
    return PL;
  }

  VEC3D_T xPower = -D1; // x~1 regime

  VEC3D_T Ln = log(D1 + x);
  VEC3D_T PL = log(x) * Ln - Ln * Ln / D2 + m_PI2d6;
  VEC3D_T dPL = m_PI2d6;
  VEC3D_T i = D0;
  VEC3D_T y = D1 / (D1 + x);

  while (fabs(dPL) > m_DE_int * m_PI2d6) {
    i += D1;
    xPower *= y;
    dPL = xPower / i / i;
    PL += dPL;
  }
  // std::cout<<" 3 "<<PL<<std::endl;
  return PL;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} // namespace IGCascade
