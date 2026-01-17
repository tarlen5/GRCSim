/*!
  \file     PairProduction.hpp
            PairProduction class header file

  \author   Tim Arlen              \n
            UCLA                   \n
      arlen@astro.ucla.edu   \n

  \author   Vladimir Vassiliev     \n
            UCLA                   \n
      vvv@astro.ucla.edu     \n

  \date     07/04/2007
  \version  0.0
  \note
*/

#ifndef IGCASCADE_PAIRPRODUCTION_H
#define IGCASCADE_PAIRPRODUCTION_H

#include "DIRBRBase.hpp"
#include "HighPrecProp.hpp"
#include "PhysicsConstants.hpp"
#include "Table.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "convert.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <qd/dd_real.h>

namespace IGCascade {

class PairProduction {
public:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors///////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Overloaded Constructor
  /// \param: RandomNumbers reference
  PairProduction(RandomNumbers *_rng, VEC3D_T ze);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Public Member Functions////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bool CheckPairProductionEBL(
      DIRBRBase *ebl_model, const double gam_ph_egy, const VEC3D_T z_o,
      const double z_min, VEC3D_T &z_int, double &TotalLambdaInt
  );
  // void GetOpticalDepth(DIRBRBase* ebl_model,const double gam_ph_egy,
  //			 const VEC3D_T z_o, VEC3D_T& z_int);
  bool UpdateGammaPhoton(
      Vec4D &gam_ph_p4, Vec4D &gam_ph_r4, VEC3D_T &gam_ph_z,
      VEC3D_T &gam_ph_z_s, VEC3D_T &delta_z_step
  );
  VEC3D_T GetEBLPhotonEgy(
      DIRBRBase *ebl_model, const double TotalLambdaInt, const double z_int,
      const double z_o, const double gam_ph_egy
  );
  void UpdateEBLPhoton(
      Vec4D &gam_ph_p4, Vec4D &gam_ph_r4, VEC3D_T &gam_ph_z,
      VEC3D_T &gam_ph_z_s, Vec4D &bg_ph_p4, Vec4D &bg_ph_r4, VEC3D_T &bg_ph_z
  );
  bool RelativisticKinematics(
      Vec4D &gam_ph_p4, Vec4D &bg_ph_p4, Vec4D &elec_p4, Vec4D &pos_p4
  );
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Public Data Members////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // computation accuracy
  VEC3D_T m_DE;

private:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Member Functions///////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Default constructor
  /// \param: N/A
  inline PairProduction();

  void DefineNumericConst(void);
  void DefineIntegrationParameters(VEC3D_T ze);

  bool IntegrateToTau(
      DIRBRBase *ebl_model, const double gam_ph_egy, const VEC3D_T z_o,
      const double z_min, const double tauFinal, VEC3D_T &z_int,
      double &TotalLambdaInt
  );

  VEC3D_T Scattering(VEC3D_T x);
  VEC3D_T PolyLog1(VEC3D_T x);

  // void NonRadialPropagation(Vec4D gam_ph_p4,Vec4D gam_ph_r4,
  //                          VEC3D_T gam_ph_z, VEC3D_T gam_ph_z_s,
  //                          VEC3D_T delta_z, VEC3D_T& L_prop,
  //                         VEC3D_T& deltaz_s, VEC3D_T& Delta_time);
  VEC3D_T ImpactAngle(VEC3D_T q);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Data Members///////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  RandomNumbers *m_rng;

  //~~~~~~~~VEC3D_T Values ~~~~~~~~~
  VEC3D_T m_dz_time_delay;
  VEC3D_T m_dz_dd;
  // VEC3D_T m_num_steps;

  double m_trans_prob_min;

  // Relevant mathematical constants
  VEC3D_T m_PI2d6;
  // VEC3D_T PI;

  // For PropagationLengthEBL integral
  // double qmin;
  double m_dlog_q;
  double m_dz;
  int m_imax;
  // double z_min;

  // double constants for increased speed
  double m_OmegaR_d;
  double m_OmegaM_d;
  double m_OmegaL_d;
  double m_Omega0_d;

  // dd_real numbers
  VEC3D_T D0;
  VEC3D_T D1;
  VEC3D_T D2;
  VEC3D_T D3;
  VEC3D_T D4;
  VEC3D_T D5;
  VEC3D_T D6;
  VEC3D_T D7;
  VEC3D_T D8;

  std::vector<double> lambda;
  std::vector<double> m_q_vec;
  std::vector<double> m_F_vec;

  // computation accuracy
  VEC3D_T m_DE_int;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Default class constructor
inline PairProduction::PairProduction() {
  // m_rng = new TRandom3(0);//("random_numbers.seed");
  m_rng = new RandomNumbers(0.0, 1.0);
}

} // namespace IGCascade
#endif // IGCASCADE_PAIRPRODUCTION_H
