/*! \file RelParticle.hpp
          RelParticle class header file

  \author   Yusef Shafi            \n
            UCLA                   \n
      yshafi@ucla.edu        \n

  \author   Tim Arlen              \n
            UCLA                   \n
      arlen@astro.ucla.edu   \n

  \date     08/07/2006

  \version  1.2

  \Revision 03/15/2008 Added redshift as a parameter

  \Revision 04/24/2008 Added electric charge as parameter

  \note  Uses typedef VEC3D_T, name changed for standardization
*/

#ifndef IGCASCADE_RELPARTICLE_H
#define IGCASCADE_RELPARTICLE_H

#ifndef __VS_NO_IOSTREAM
#include <iostream>
#endif

#include <cmath>
#include <string>
#include <vector>

#include <qd/dd_real.h>

#include "Vec3D.hpp"
#include "Vec4D.hpp"

namespace IGCascade {
class RelParticle {

public:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors///////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Default constructor
  inline RelParticle();

  // Overloaded Constructors
  inline RelParticle(Vec4D &_r4, Vec4D &_p4);
  inline RelParticle(Vec4D &_r4, Vec4D &_p4, VEC3D_T _z);
  inline RelParticle(Vec4D &_r4, Vec4D &_p4, VEC3D_T _z, int _q);
  inline RelParticle(VEC3D_T _p0, VEC3D_T _m0, Vec4D &_r4, Vec3D &_n);

  inline RelParticle(RelParticle &_rp);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Static Functions///////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Find/Boost to Center of Momentum Frame:
  static Vec4D GetUcm(std::vector<RelParticle> &particles);
  static void BoostSystem(Vec4D &u4, std::vector<RelParticle> &particles);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Member Functions///////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // When the particle is boosted, just call Vec4D's two boosts
  inline void Boost(Vec4D &u4);     // particle boost
  inline void Rotate(Vec3D &axis);  // particle rotation
  inline void Reflect(Vec3D &norm); // particle reflection
  inline void P();                  // parity transformation
  inline void T();                  // time inversion

  // Dump Function to get info out
  // void Dump(std::ostream& stream = std::cout) const;

  // void Propagate() - give a dt and velocity using += operators, etc
  inline void PropagateFree(VEC3D_T dt);

  // Initialize all data members to zero...
  inline void Zero();

  // choose a random vector on sphere (uniform distribution)
  /***Moved to Vec3D***/
  /*Vec3D Scatter3DUniform();*/

  void Dump(std::ostream &stream = std::cout) const; //!< prints coordinates
  void DumpShort(std::ostream &stream = std::cout) const;
  //!< prints coordinates

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Public Data Members////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Vec4D m_r4;          // Position
  Vec4D m_p4;          // Energy-Momentum
  VEC3D_T m_m0;        // Mass
  VEC3D_T m_z;         // Redshift of particle
  VEC3D_T m_z_s;       // Redshift of test particle at current distance
  VEC3D_T m_egy_int;   // Energy of Photon at IC creation point
  VEC3D_T m_z_s_int;   // z_s of GammaPhoton at IC point (when created)
  VEC3D_T m_z_int;     // z of Gamma at IC point (when created)
  VEC3D_T m_rnext;     // Radial coordinate (next)
  VEC3D_T m_elec_time; // another temporary needed variable, to
                       // track time delay due exclusively to
                       // electrons
  int m_q; // Charge in units of |e| (i.e. electron has q = -1)
  int m_tag;
  double m_weight; // weight of each particle (bc small tr prob correction)

  //////////////////////////////////////////////////////////
  //////////// 9 Branches for Track Lepton Tree ////////////
  double m_tz;
  double m_tegy;
  double m_tpx;
  double m_tpy;
  double m_tpz;
  double m_ttime;
  double m_trx;
  double m_try;
  double m_trz;

private:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Data Members///////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // std::string m_tag;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Default class constructor
inline RelParticle::RelParticle() : m_r4(), m_p4() {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor (1)
/// \param _r: r
/// \param _ct: ct
inline RelParticle::RelParticle(Vec4D &_r4, Vec4D &_p4) {
  m_r4 = _r4;
  m_p4 = _p4;
  m_m0 = sqrt(m_p4.Norm2());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor (2)
/// \param _r4: 4 position
/// \param _p4: 4 momentum
/// \param _z:  redshift
inline RelParticle::RelParticle(Vec4D &_r4, Vec4D &_p4, VEC3D_T _z) {
  m_r4 = _r4;
  m_p4 = _p4;
  m_z = _z;
  m_m0 = sqrt(m_p4.Norm2());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor (3)
/// \param _r4: 4 position
/// \param _p4: 4 momentum
/// \param _z:  redshift
/// \param _q:  charge
inline RelParticle::RelParticle(Vec4D &_r4, Vec4D &_p4, VEC3D_T _z, int _q) {
  m_r4 = _r4;
  m_p4 = _p4;
  m_z = _z;
  // If the mass is supposed to be zero, it may turn negative, because
  //   we lose precision on the computation for Norm2().
  if (m_p4.Norm2() < 0.0)
    m_m0 = 0.0;
  else
    m_m0 = sqrt(m_p4.Norm2());
  m_q = _q;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor (4)
/// \param _p0: Energy  [eV]
/// \param _m0: rest mass [eV]
/// \param _n: direction
inline RelParticle::RelParticle(VEC3D_T _p0, VEC3D_T _m0, Vec4D &_r4,
                                Vec3D &_n) {
  m_r4 = _r4;
  m_m0 = _m0;
  // m_p4 = Vec4D(_p0, sqrt(_p0*_p0 - m_m0*m_m0)*_n);
  m_p4.r0 = _p0;
  m_p4.r = _n;
  m_p4.r *= sqrt(_p0 * _p0 - m_m0 * m_m0);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Copy constructor
// \param _rp: RelParticle object to copy
inline RelParticle::RelParticle(RelParticle &_rp)
    : m_r4(_rp.m_r4), m_p4(_rp.m_p4) {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Member functions ///
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Boost function
// \param: 4D velocity
inline void RelParticle::Boost(Vec4D &u4) {
  m_r4.Boost4D(u4);
  m_p4.Boost4D(u4);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Rotate function
// \param: axis about which to rotate
inline void RelParticle::Rotate(Vec3D &axis) {
  m_r4.Rotate(axis);
  m_p4.Rotate(axis);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reflect function
// \param: normalized axis about which to rotate
inline void RelParticle::Reflect(Vec3D &axis) {
  m_r4.Reflect(axis);
  m_p4.Reflect(axis);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Parity Transformation function
// \param: N/A
inline void RelParticle::P() {
  m_r4.P();
  m_p4.P();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Time inversion function
// \param: N/A
inline void RelParticle::T() {
  m_r4.T();
  m_p4.T();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialization to zero
// \param: N/A
inline void RelParticle::Zero() {
  Vec4D zero_vec(0., 0., 0., 0.);
  m_r4 = zero_vec;
  m_p4 = zero_vec;
  m_m0 = 0.0;
  m_z = 0.0;
  m_z_s = 0.0;
  m_z_s_int = 0.0;
  m_z_int = 0.0;
  m_rnext = 0.0;
  m_q = 0;
  m_tag = 0;
  m_weight = 0.0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Propagate function
// \param: time increment
inline void RelParticle::PropagateFree(VEC3D_T dt) {

  m_r4.r.x += (m_p4.r.x / m_p4.r0) * dt;
  m_r4.r.y += (m_p4.r.y / m_p4.r0) * dt;
  m_r4.r.z += (m_p4.r.z / m_p4.r0) * dt;
  m_r4.r0 += dt;
}

} // namespace IGCascade

#ifndef __VS_NO_IOSTREAM

namespace std {
inline ostream &operator<<(ostream &stream, IGCascade::RelParticle &R);
inline istream &operator>>(istream &stream, IGCascade::RelParticle &R);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream insertion
inline ostream &operator<<(ostream &stream, IGCascade::RelParticle &R) {
  // R.DumpShort(stream);
  R.Dump(stream);
  return stream;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream extraction
inline istream &operator>>(istream &stream, IGCascade::RelParticle &R) {
  char c;
  stream >> c >> R.m_r4.r.x >> R.m_r4.r.y >> R.m_r4.r.z >> c;
  return stream;
}
} // namespace std

#endif

#endif // IGCASCADE_RELPARTICLE_H
