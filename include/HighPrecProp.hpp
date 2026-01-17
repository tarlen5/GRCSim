/*!  HighPrecProp.hpp

  author: Timothy C. Arlen
          tca3@psu.edu

  date:    24 October 2014

  A set of functions that do the high-precision propagation, updating
  the particle's (electron or photon) kinematic parameters throuhgout
  Cosmologically expanding spacetime.

 */

#ifndef HIGHPRECPROP_H
#define HIGHPRECPROP_H

#include <cmath>

#include <qd/dd_real.h>

#include "PhysicsConstants.hpp"
#include "RelParticle.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"

using IGCascade::RelParticle;
using IGCascade::Vec3D;
using IGCascade::Vec4D;

using PhysConst::CGS_C;
using PhysConst::CGS_HUBRAD;
using PhysConst::HUB_CONST;
using PhysConst::OMEGA_0;
using PhysConst::OMEGA_L;
using PhysConst::OMEGA_M;
using PhysConst::OMEGA_R;

inline void NonRadialPhotonPropagation(Vec4D gam_ph_p4, Vec4D gam_ph_r4,
                                       VEC3D_T gam_ph_z, VEC3D_T gam_ph_z_s,
                                       VEC3D_T delta_z, VEC3D_T &L_prop,
                                       VEC3D_T &deltaz_s, VEC3D_T &Delta_time)
/*!
  For a non-radially propagating photon, calculates updated dynamical
  parameters to fourth order accuracy.

  \param
      gam_ph_p4  - gamma ray 4 momentum
      gam_ph_r4  - gamma ray 4 position
      gam_ph_z   - gamma ray redshift
      gam_ph_z_s - gamma ray z_s, ( ~ time delay; see Cosmology writeup!)
      L_prop     - dist traveled by gam_ph; calculated here
      delta_z    - used in calculations (NOT updated here)
      deltaz_s   - calculated here
      Delta_time - calculated here to fourth order.

  \return - none

*/
{
  VEC3D_T z_dd = gam_ph_z;
  VEC3D_T z_s = gam_ph_z_s;

  VEC3D_T delta_z2 = delta_z * delta_z;
  VEC3D_T delta_z3 = delta_z * delta_z2;
  VEC3D_T delta_z4 = delta_z * delta_z3;

  VEC3D_T z1 = (1. + z_dd);
  VEC3D_T z2 = z1 * z1;
  VEC3D_T z3 = z1 * z2;
  VEC3D_T z4 = z1 * z3;
  VEC3D_T Q =
      sqrt(OMEGA_R * z4 + OMEGA_M * z3 + OMEGA_L + (1.0 - OMEGA_0) * z2);
  VEC3D_T a = 1.0 / Q;
  VEC3D_T a2 = a * a;
  VEC3D_T a3 = a * a2;
  VEC3D_T a4 = a * a3;
  VEC3D_T G =
      2.0 * OMEGA_R * z3 + 3.0 / 2.0 * OMEGA_M * z2 + (1.0 - OMEGA_0) * z1;
  VEC3D_T H = 6.0 * OMEGA_R * z2 + 3.0 * OMEGA_M * z1 + (1.0 - OMEGA_0);
  VEC3D_T J = 12.0 * OMEGA_R * z1 + 3.0 * OMEGA_M;
  VEC3D_T b = -a3 * G;
  VEC3D_T b2 = b * b;
  VEC3D_T b3 = b * b2;
  VEC3D_T c = -3.0 * a2 * b * G - a3 * H;
  VEC3D_T d = -3.0 * a * (2.0 * b2 + a * c) * G - 6.0 * a2 * b * H - a3 * J;

  L_prop = a * delta_z - b * delta_z2 / 2.0 + c * delta_z3 / 6.0 -
           d * delta_z4 / 24.0;
  VEC3D_T L_prop2 = L_prop * L_prop;
  VEC3D_T L_prop3 = L_prop * L_prop2;
  VEC3D_T L_prop4 = L_prop * L_prop3;

  VEC3D_T alpha = a / z1;
  VEC3D_T beta = -a / z2 + b / z1;
  VEC3D_T gamma = 2.0 * a / z3 - 2.0 * b / z2 + c / z1;
  VEC3D_T delta = -6.0 * a / z4 + 6.0 * b / z3 - 3.0 * c / z2 + d / z1;

  VEC3D_T Tau = alpha / a * L_prop +
                (b * alpha / a - beta) / a2 / 2.0 * L_prop2 +
                (gamma / 6.0 - b * beta / 2.0 / a - c * alpha / 6.0 / a +
                 b2 * alpha / 2.0 / a2) /
                    a3 * L_prop3 +
                (b * gamma / 4.0 / a - delta / 24.0 + c * beta / 6.0 / a +
                 d * alpha / 24.0 / a - 5.0 / 12.0 * b * c * alpha / a2 -
                 5.0 / 8.0 * b2 * beta / a2 + 5.0 / 8.0 * b3 * alpha / a3) /
                    a4 * L_prop4;

  Vec3D reH = gam_ph_r4.r / PhysConst::CGS_HUBRAD;
  Vec3D e_L = gam_ph_p4.r / gam_ph_p4.r.Norm();
  VEC3D_T L_prop_s =
      L_prop * (2.0 * reH * e_L + L_prop) /
      (sqrt(reH * reH + 2.0 * L_prop * reH * e_L + L_prop * L_prop) +
       reH.Norm());
  VEC3D_T L_prop_s2 = L_prop_s * L_prop_s;
  VEC3D_T L_prop_s3 = L_prop_s * L_prop_s2;
  VEC3D_T L_prop_s4 = L_prop_s * L_prop_s3;

  z1 = (1. + z_s);
  z2 = z1 * z1;
  z3 = z1 * z2;
  z4 = z1 * z3;
  Q = sqrt(OMEGA_R * z4 + OMEGA_M * z3 + OMEGA_L + (1.0 - OMEGA_0) * z2);
  a = 1.0 / Q;
  a2 = a * a;
  a3 = a * a2;
  a4 = a * a3;
  G = 2.0 * OMEGA_R * z3 + 3.0 / 2.0 * OMEGA_M * z2 + (1.0 - OMEGA_0) * z1;
  H = 6.0 * OMEGA_R * z2 + 3.0 * OMEGA_M * z1 + (1.0 - OMEGA_0);
  J = 12.0 * OMEGA_R * z1 + 3.0 * OMEGA_M;
  b = -a3 * G;
  b2 = b * b;
  b3 = b * b2;
  c = -3.0 * a2 * b * G - a3 * H;
  d = -3.0 * a * (2.0 * b2 + a * c) * G - 6.0 * a2 * b * H - a3 * J;

  alpha = a / z1;
  beta = -a / z2 + b / z1;
  gamma = 2.0 * a / z3 - 2.0 * b / z2 + c / z1;
  delta = -6.0 * a / z4 + 6.0 * b / z3 - 3.0 * c / z2 + d / z1;

  VEC3D_T Tau_s = alpha / a * L_prop_s +
                  (b * alpha / a - beta) / a2 / 2.0 * L_prop_s2 +
                  (gamma / 6.0 - b * beta / 2.0 / a - c * alpha / 6.0 / a +
                   b2 * alpha / 2.0 / a2) /
                      a3 * L_prop_s3 +
                  (b * gamma / 4.0 / a - delta / 24.0 + c * beta / 6.0 / a +
                   d * alpha / 24.0 / a - 5.0 / 12.0 * b * c * alpha / a2 -
                   5.0 / 8.0 * b2 * beta / a2 + 5.0 / 8.0 * b3 * alpha / a3) /
                      a4 * L_prop_s4;

  Delta_time = (Tau - Tau_s) / PhysConst::HUB_CONST;

  VEC3D_T A = 1.0 / a;
  VEC3D_T B = b / 2.0 / a3;
  VEC3D_T C = (b2 / 2.0 / a2 - c / 6.0 / a) / a3;
  VEC3D_T D =
      (d / 24.0 - 5.0 / 12.0 * b * c / a + 5.0 / 8.0 * b3 / a2) / a3 / a2;

  deltaz_s = A * L_prop_s + B * L_prop_s2 + C * L_prop_s3 + D * L_prop_s4;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// FROM GALACTIC GRID:

inline VEC3D_T GetLeptonDelta_z(VEC3D_T &PL, RelParticle *&Lepton)
/*!
  Calculates change in redshift for a lepton propagating a distance PL
  IMPORTANT: does not change any of Lepton's values.

  \param PL     - prop length [cm]
  \param Lepton - charged lepton which is created by IC process

  \return delta_z - change in redshift of lepton

  **NOTE** Changed: E -> E_i *( 1+z_{i+1} )/( 1+z_i) in redshift expansion
  */
{

  VEC3D_T m_eV = sqrt(Lepton->m_p4 * Lepton->m_p4);
  VEC3D_T ze = Lepton->m_z;
  VEC3D_T p_o = Lepton->m_p4.r.Norm() / (1.0 + ze);
  VEC3D_T rho = m_eV * m_eV / p_o / p_o;

  VEC3D_T z1 = (1.0 + ze);
  VEC3D_T z2 = z1 * z1;
  VEC3D_T z3 = z1 * z2;
  VEC3D_T z4 = z1 * z3;
  VEC3D_T Q =
      sqrt(OMEGA_R * z4 + OMEGA_M * z3 + OMEGA_L + (1.0 - OMEGA_0) * z2);

  VEC3D_T F = 1.0 / Q;
  VEC3D_T F2 = F * F;
  VEC3D_T F3 = F * F2;
  VEC3D_T G =
      2.0 * OMEGA_R * z3 + 3.0 / 2.0 * OMEGA_M * z2 + (1.0 - OMEGA_0) * z1;
  VEC3D_T H = 6.0 * OMEGA_R * z2 + 3.0 * OMEGA_M * z1 + (1.0 - OMEGA_0);
  // VEC3D_T J  = 12.0*OMEGA_R*z1 + 3.0*OMEGA_M;
  VEC3D_T F_pr1 = -F3 * G;
  VEC3D_T F_pr2 = -3.0 * F2 * F_pr1 * G - F3 * H;
  // VEC3D_T F_pr3= -3.0*F*(2.0*F_pr1*F_pr1 + F*F_pr2)*G - 6.0*F2*F_pr1*H -F3*J;

  VEC3D_T P = (1.0 - rho / z2 / 2.0);
  VEC3D_T P_pr1 = rho / z3;
  VEC3D_T P_pr2 = -3.0 * rho / z4;
  // VEC3D_T P_pr3 = 12.0*rho/(z2*z3);

  VEC3D_T a = P * F;
  VEC3D_T a2 = a * a;
  VEC3D_T a3 = a * a2;
  // VEC3D_T a4 = a*a3;
  VEC3D_T b = P_pr1 * F + P * F_pr1;
  VEC3D_T b2 = b * b;
  // VEC3D_T b3 = b*b2;
  VEC3D_T c = P_pr2 * F + 2.0 * P_pr1 * F_pr1 + P * F_pr2;
  // VEC3D_T d = P_pr3*F + 3.0*P_pr2*F_pr1 + 3.0*P_pr1*F_pr2 + P*F_pr3;

  VEC3D_T A = 1.0 / a;
  VEC3D_T B = b / 2.0 / a3;
  VEC3D_T C = (b2 / 2.0 / a2 - c / 6.0 / a) / a3;
  // VEC3D_T D = (d/24.0 - 5.0/12.0*b*c/a + 5.0/8.0*b3/a2)/a3/a2;

  VEC3D_T L_prop = PL / CGS_HUBRAD;
  VEC3D_T L_prop2 = L_prop * L_prop;
  VEC3D_T L_prop3 = L_prop * L_prop2;
  // VEC3D_T L_prop4 = L_prop*L_prop3;

  VEC3D_T delta_z = A * L_prop + B * L_prop2 + C * L_prop3; // + D*L_prop4;

  return delta_z;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline VEC3D_T GetPhotonDelta_zs(VEC3D_T &PL, RelParticle *&Lepton,
                                 Vec3D &r_new, VEC3D_T &delta_time)
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

  VEC3D_T m_eV = sqrt(Lepton->m_p4 * Lepton->m_p4);
  // std::cerr<<"m_eV (GetLeptonGetPhotonDelta_zs): "<<m_eV<<std::endl;
  VEC3D_T ze = Lepton->m_z;
  VEC3D_T p_o = Lepton->m_p4.r.Norm() / (1.0 + ze);
  VEC3D_T rho = m_eV * m_eV / p_o / p_o;

  //------------------------ BEGIN z_s Computation ------------------------//
  VEC3D_T z = Lepton->m_z_s;
  VEC3D_T z1 = (1. + z);
  VEC3D_T z2 = z1 * z1;
  VEC3D_T z3 = z2 * z1;
  VEC3D_T z4 = z3 * z1;
  VEC3D_T Q =
      sqrt(OMEGA_R * z4 + OMEGA_M * z3 + OMEGA_L + (1.0 - OMEGA_0) * z2);

  VEC3D_T a = 1.0 / Q;
  VEC3D_T a2 = a * a;
  VEC3D_T a3 = a * a2;
  // VEC3D_T a4 = a*a3;
  VEC3D_T G =
      2.0 * OMEGA_R * z3 + 3.0 / 2.0 * OMEGA_M * z2 + (1.0 - OMEGA_0) * z1;
  VEC3D_T H = 6.0 * OMEGA_R * z2 + 3.0 * OMEGA_M * z1 + (1.0 - OMEGA_0);
  VEC3D_T J = 12.0 * OMEGA_R * z1 + 3.0 * OMEGA_M;
  VEC3D_T b = -a3 * G;
  VEC3D_T b2 = b * b;
  // VEC3D_T b3 = b*b2;
  VEC3D_T c = -3.0 * a2 * b * G - a3 * H;
  // VEC3D_T d  = -3.0*a*(2.0*b2 + a*c)*G - 6.0*a2*b*H - a3*J;

  VEC3D_T A = 1.0 / a;
  VEC3D_T B = b / 2.0 / a3;
  VEC3D_T C = (b2 / 2.0 / a2 - c / 6.0 / a) / a3;
  // VEC3D_T D = (d/24.0 - 5.0/12.0*b*c/a + 5.0/8.0*b3/a2)/a3/a2;

  VEC3D_T L_prop_s = (r_new.Norm() - Lepton->m_r4.r.Norm()) / CGS_HUBRAD;
  VEC3D_T L_prop_s2 = L_prop_s * L_prop_s;
  VEC3D_T L_prop_s3 = L_prop_s * L_prop_s2;
  // VEC3D_T L_prop_s4 = L_prop_s*L_prop_s3;

  VEC3D_T delta_z_s =
      A * L_prop_s + B * L_prop_s2 + C * L_prop_s3; //+ D*L_prop_s4;
  VEC3D_T delta_z_s2 = delta_z_s * delta_z_s;
  VEC3D_T delta_z_s3 = delta_z_s * delta_z_s2;
  // VEC3D_T delta_z_s4 = delta_z_s*delta_z_s3;

  VEC3D_T alpha = a / z1;
  VEC3D_T beta = -a / z2 + b / z1;
  VEC3D_T gamma = 2.0 * a / z3 - 2.0 * b / z2 + c / z1;
  // VEC3D_T delta = -6.0*a/z4 + 6.0*b/z3 - 3.0*c/z2 + d/z1;

  VEC3D_T tau_s = alpha * delta_z_s - beta / 2.0 * delta_z_s2 +
                  gamma / 6.0 * delta_z_s3; // - delta/24.0*delta_z_s4;

  // Time Delay Testing:
  // std::ofstream TimeList("time_cmds.txt",std::ios::app);

  //     std::string cmd_begin = "tau_strue =
  //     integrate[1/(1+z)/sqrt[OmegaR*((1+z)^4) + OmegaM*((1+z)^3) + OmegaL +
  //     (1-Omega0)*((1+z)^2)], \\ \n";
  //     TimeList<<std::setprecision(32)<<cmd_begin<<"{z, "
  //          <<(z-delta_z_s)<<", "<<z<< "}]"<<std::endl;
  //     TimeList<<"tau_s = "<<tau_s<<std::endl;
  //     TimeList<<"(tau_s - tau_strue)/tau_strue"<<std::endl;

  //     std::cout<<(z - delta_z_s)<<" "<<z<<std::endl;

  //     TimeList.close();

  //     char getline;
  //     std::cin>>getline;

  //------------------------ END z_s -----------------------------//

  //------------------------ BEGIN z Computation --------------------------//
  VEC3D_T delta_z = GetLeptonDelta_z(PL, Lepton);
  VEC3D_T delta_z2 = delta_z * delta_z;
  VEC3D_T delta_z3 = delta_z * delta_z2;
  // VEC3D_T delta_z4 = delta_z*delta_z3;

  // Compute Phi Derivatives:
  z = Lepton->m_z;
  z1 = (1. + z);
  z2 = z1 * z1;
  z3 = z2 * z1;
  z4 = z3 * z1;
  Q = sqrt(OMEGA_R * z4 + OMEGA_M * z3 + OMEGA_L + (1.0 - OMEGA_0) * z2);
  VEC3D_T F = 1.0 / Q;
  VEC3D_T F2 = F * F;
  VEC3D_T F3 = F * F2;
  G = 2.0 * OMEGA_R * z3 + 3.0 / 2.0 * OMEGA_M * z2 + (1.0 - OMEGA_0) * z1;
  H = 6.0 * OMEGA_R * z2 + 3.0 * OMEGA_M * z1 + (1.0 - OMEGA_0);
  J = 12.0 * OMEGA_R * z1 + 3.0 * OMEGA_M;
  VEC3D_T F_pr1 = -F3 * G;
  VEC3D_T F_pr2 = 3.0 / 2.0 * G * F3 * F2 - F3 * H;
  // VEC3D_T F_pr3= -3.0*F*(2.0*F_pr1*F_pr1 + F*F_pr2)*G-6.0*F2*F_pr1*H-F3*J;
  // VEC3D_T F_pr3 = F3*(3.0/2.0*(H+J)*F2 - 15.0/4.0*G*F3*F - J );

  VEC3D_T P = (1.0 - rho / z2 / 2.0);
  VEC3D_T K = P / z1;
  VEC3D_T K_pr = (3.0 / 2.0 * rho / z2 - 1.0) / z2;
  VEC3D_T K_pr2 = (2.0 - 6.0 * rho / z2) / z3;
  // VEC3D_T K_pr3 = 6.0*(5.0*rho/z2 - 1.0)/z4;

  VEC3D_T Phi = K * F;
  VEC3D_T Phi_pr1 = K_pr * F + K * F_pr1;
  VEC3D_T Phi_pr2 = K_pr2 * F + 2.0 * K_pr * F_pr1 + K * F_pr2;
  // VEC3D_T Phi_pr3 = K_pr3*F + 3.0*K_pr2*F_pr1 + 3.0*K_pr*F_pr2 + K*F_pr3;

  VEC3D_T tau = Phi * delta_z - Phi_pr1 / 2.0 * delta_z2 +
                Phi_pr2 / 6.0 * delta_z3; // - Phi_pr3/24.0*delta_z4;

  // Time Delay Testing:
  //     std::ofstream TimeList("time_cmds.txt",std::ios::app);

  //     TimeList<<std::setprecision(32)<<"rho = "<<rho<<std::endl;

  //     std::string cmd_begin = "tau_true = integrate[1/sqrt[OmegaR*((1+z)^4) +
  //     OmegaM*((1+z)^3) + OmegaL + (1-Omega0)*((1+z)^2)]/sqrt[rho+(1+z)^2],
  //     \\ \n";
  //     //std::string cmd_begin = "tau_true = integrate[1/sqrt[OmegaR*((1+z)^4)
  //     + OmegaM*((1+z)^3) + OmegaL +
  //     (1-Omega0)*((1+z)^2)]/(1+z)/(1-rho/2/(1+z)^2), \\ \n";
  //     TimeList<<std::setprecision(32)<<cmd_begin<<"{z, "
  //          <<(z-delta_z)<<", "<<z<< "}]"<<std::endl;
  //     TimeList<<"tau = "<<tau<<std::endl;
  //     TimeList<<"(tau - tau_true)/tau_true"<<std::endl;

  //     std::cout<<(z - delta_z)<<" "<<z<<" "<<delta_z<<std::endl;

  //     TimeList.close();

  //     char getline;
  //     std::cin>>getline;

  //------------------------ END z -----------------------------//

  VEC3D_T time_delay = (tau - tau_s) / HUB_CONST;

  delta_time = time_delay;
  return delta_z_s;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline void PropagateConstantMF(RelParticle &Photon, RelParticle *&Lepton,
                                const VEC3D_T &mf_magnitude, VEC3D_T &PL,
                                Vec3D &r_new, VEC3D_T &time_delay, Vec3D &n_eo,
                                VEC3D_T &delta_z, VEC3D_T &delta_zs, Vec3D &e_b)
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

  VEC3D_T ze = Lepton->m_z;
  VEC3D_T m_eV = sqrt(Lepton->m_p4 * Lepton->m_p4);
  VEC3D_T p_o = Lepton->m_p4.r.Norm() / (1.0 + ze);

  VEC3D_T Kappa_o = CGS_C * 1.E-8;
  Kappa_o *= Lepton->m_q * mf_magnitude / p_o;

  VEC3D_T KappaL = Kappa_o * PL;

  // Update Lepton parameters:
  Lepton->m_r4.r = r_new;
  Lepton->m_z -= delta_z;
  Lepton->m_z_s -= delta_zs;
  Lepton->m_r4.r0 += time_delay;

  //-----------------------------------------------------------------
  // Photon rotation
  //-----------------------------------------------------------------
  Vec3D n_eo_ = (n_eo - (n_eo * e_b) * e_b) * cos(KappaL) +
                (n_eo ^ e_b) * sin(KappaL) + (n_eo * e_b) * e_b;

  // Unless rotation angle is very small, do rotation:
  Vec3D n_cp = n_eo_ ^ n_eo;

  // If theta is close to Pi/2, this will avoid the domain error if arg
  //   to asin() ~ 1.0 + epsilon.
  VEC3D_T theta = "0.0";
  if (n_cp.Norm() < 0.99)
    theta = asin(n_cp.Norm());
  else {
    VEC3D_T n_dp = (n_eo_ * n_eo) / n_eo_.Norm() / n_eo.Norm();
    theta = acos(n_dp);
  }

  //   VEC3D_T theta = asin(n_cp.Norm());
  if (fabs(theta) > 1.0E-15) {
    Vec3D axis = theta * n_cp / n_cp.Norm();
    Vec3D n_photon(Photon.m_p4.r / Photon.m_p4.r.Norm());
    n_photon.Rotate(axis);
    Photon.m_p4.r0 *= (1 + ze - delta_z) / (1 + ze);
    Photon.m_p4.r = Photon.m_p4.r0 * n_photon;
  } else // No rotation, just modify 4-mom due to expansion
  {
    Photon.m_p4.r0 *= (1 + ze - delta_z) / (1 + ze);
    Photon.m_p4.r *= Photon.m_p4.r0 / Photon.m_p4.r0;
  }

  Lepton->m_p4.r = p_o * (1.0 + Lepton->m_z) * n_eo_;
  Lepton->m_p4.r0 = sqrt(Lepton->m_p4.r * Lepton->m_p4.r + m_eV * m_eV);
  // std::cerr<<"Lepton 4 mom: "<<Lepton->m_p4<<std::endl<<std::endl;

  Photon.m_r4 = Lepton->m_r4;
  Photon.m_z = Lepton->m_z;
  Photon.m_z_s = Lepton->m_z_s;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif
