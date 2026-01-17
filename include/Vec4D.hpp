//-*-mode:c++; mode:font-lock;-*-

/*! \file Vec4D.hpp
  Vec4D class header file

  \author   Vladimir Vassiliev        \n
            UCLA                      \n
      vvv@astro.ucla.edu        \n

  \author   Yusef Shafi               \n
            UCLA                      \n
      yshafi@ucla.edu           \n

  \date     08/07/2006

  \version  1.3

  \note	    1.2 - We change all doubles to double-double precision
            using dd lib. We also add a new boost function that
      takes a velocity four vector.
      1.3 - We add a function to return a normalized direction vector.
*/

#ifndef IGCASCADE_VEC4D_H
#define IGCASCADE_VEC4D_H

#include <qd/dd_real.h>
// #include<qd/dd.h>

#include "Vec3D.hpp"

namespace IGCascade {

/*!  \class Vec4D
  \brief 4 dimensional Lorentz vector class

  This class defines 4D vectors, set of operations
   in vector field, scalar product, rotations and boosts
   as well as parity and time transformations.

*/

class Vec4D {
public:
  inline Vec4D();                      //!< default constructor
  inline Vec4D(const Vec4D &v);        //!< copy constructor
  inline Vec4D(VEC3D_T _r0, Vec3D _r); //!< constructor
  inline Vec4D(VEC3D_T _r0, VEC3D_T r1, VEC3D_T r2,
               VEC3D_T r3); //!< constructor

  inline VEC3D_T Norm2() const; //!< calculates scalar product with itself

  // Original transformation operations
  inline void Rotate(const Vec3D &axis); //!< rotates vector around axis
  inline bool Boost(const Vec3D &boost); //!< boosts vector
  inline Vec3D SRFBoost() const; //!< finds boost to special reference frame
  inline void P();               //!< parity transformation
  inline void T();               //!< time inversion
  inline void Reflect(const Vec3D &norm); //!< reflect vector

  // Additional transformation operations
  inline bool Boost4D(const Vec4D &boost4D); //!< boosts vector with Vec4D

  // Obtain normalized direction vector
  inline Vec3D GetDirection();

  inline Vec4D &Reset(const Vec4D &v = Vec4D());
  inline Vec4D &Reset(VEC3D_T _r0, VEC3D_T r1, VEC3D_T r2, VEC3D_T r3);

  inline Vec4D &operator=(const Vec4D &v);  //!< assignment
  inline Vec4D &operator+=(const Vec4D &v); //!< assignment: addition
  inline Vec4D &operator-=(const Vec4D &v); //!< assignment: subtraction
  inline Vec4D &operator*=(VEC3D_T d);      //!< assignment: multiply by scaler
  inline Vec4D &operator/=(VEC3D_T d);      //!< assignment: divide by scaler

  inline Vec4D operator+(const Vec4D &v) const;   //!< addition
  inline Vec4D operator-(const Vec4D &v) const;   //!< subtraction
  inline VEC3D_T operator*(const Vec4D &v) const; //!< scalar product

  void Dump(std::ostream &stream = std::cout) const; //!< prints coordinates
  void DumpShort(std::ostream &stream = std::cout) const; //!< prints
                                                          //!< coordinates

public:
  VEC3D_T r0; //!< time component
  Vec3D r;    //!< space components
};

inline Vec4D operator-(const Vec4D &v);            //!< negation
inline Vec4D operator*(VEC3D_T d, const Vec4D &v); //!< left scalar mult.
inline Vec4D operator*(const Vec4D &v, VEC3D_T d); //!< right scalar mult.
inline Vec4D operator/(const Vec4D &v, VEC3D_T d); //!< right scalar division

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Default class constructor
inline Vec4D::Vec4D() : r0(), r() {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Copy constructor
inline Vec4D::Vec4D(const Vec4D &v) : r0(v.r0), r(v.r) {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _r0: time component
/// \param _r: space components
inline Vec4D::Vec4D(VEC3D_T _r0, Vec3D _r) : r0(_r0), r(_r) {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor
/// \param _r0: time component
/// \param r1: space components in x-direction
/// \param r2: space components in y-direction
/// \param r3: space components in z-direction

inline Vec4D::Vec4D(VEC3D_T _r0, VEC3D_T r1, VEC3D_T r2, VEC3D_T r3)
    : r0(_r0), r(r1, r2, r3) {
  // nothing to see here
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Member functions ///
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector norm in square
inline VEC3D_T Vec4D::Norm2() const { return r0 * r0 - r * r; }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector rotation around axis vector.
/// \param axis: axis of rotation vector [rad]
/*! \note
  The modulus of rotation angle is equal to the axis
  norm and is given in radians. The rotation of e_x
  around e_z by PI/2 is equal to e_y.
*/
inline void Vec4D::Rotate(const Vec3D &axis) {
  r.Rotate(axis);
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector boost
/// \param boost: boost vector (boost,boost) < 1.0
/*! \note

*/

// Original Boost
inline bool Vec4D::Boost(const Vec3D &boost) {
  VEC3D_T s = boost.Norm();

  if (s >= (VEC3D_T)1.) {
    return false; // boost is invalid
  }

  Vec3D e_boost;
  VEC3D_T gamma = 1. / sqrt(1. - (boost * boost));
  if (s != (VEC3D_T)0) {
    e_boost = boost / s;
  } else {
    return true; // no boost required, boost=0..
  }

  VEC3D_T par = (r * e_boost);
  Vec4D t(gamma * (r0 - (boost * r)),
          gamma * (par * e_boost - r0 * boost) + r - par * e_boost);

  *this = t;

  if (gamma >= 1.0E+6) {
    std::cout << "Boost: Boost accuracy may be lost " << std::endl;
    std::cout << "       Boost is too large > 1E+6 " << std::endl;
  }

  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector boost using Vec4D
/// \param boost4D: boost4D Vec4D boost,gamma (boost,boost) < 1.0
/*! \note

*/

// new Vec4D that takes Vec4D
// note no further computation of gamma is required
inline bool Vec4D::Boost4D(const Vec4D &boost4D) {
  /*Here we begin by taking Vec4D and extracting
  beta-gamma boost vector parameter as Vec3D
  boost after dividing by gamma and last
  coordinate as scalar gamma*/

  // VEC3D_T _n = Norm2();

  VEC3D_T gamma = boost4D.r0;
  Vec3D boost = (boost4D.r) / gamma;

  VEC3D_T s = boost.Norm();

  if (s > (VEC3D_T)1) {
    return false;
  }

  Vec3D e_boost;

  if (s != (VEC3D_T)0) {
    e_boost = boost / s;
  } else {
    return true; // no boost required, boost=0..
  }

  // change the boosting specifics by factoring to preserve precision
  VEC3D_T par = (r * e_boost);
  Vec4D t(gamma * r0 - (boost4D.r) * r,
          r - par * e_boost + gamma * (par * e_boost) - (boost4D.r) * r0);

  /*VEC3D_T _n2 = t.Norm2();
  std::cout << _n2;
  if(_n>_n2)
  {
    t.r0+=sqrt(_n-_n2);
  }
  else if(_n<_n2)
  {
    t.r0-=sqrt(_n-_n2);
  }*/

  *this = t;

  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method to return a normalized direction vector
inline Vec3D Vec4D::GetDirection() {
  if (r.Norm2() == 0.) {
    return r;
  } else {
    VEC3D_T norm = r.Norm();
    Vec3D r1 = r / norm;
    return r1;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for finding boost to special r.f.
inline Vec3D Vec4D::SRFBoost() const {
  VEC3D_T n = Norm2();
  if (n > 0.) { // time-like vector
    // e.g. physical particle

    return r / r0; // boost to r.f. where r=0
                   // i.e. particle's rest frame

  } else if (n < 0.) { // space-like vector
    // e.g. tachyon

    return r0 * r / (r * r); // boost to r.f. where r0=0
  }
  // light-like
  return Vec3D(); // boost=0. is returned
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for parity transformation of Lorentz vector
inline void Vec4D::P() {
  r = -r;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for time inversion of Lorentz vector
inline void Vec4D::T() {
  r0 = -r0;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for reflection in mirror (defined by norm)
inline void Vec4D::Reflect(const Vec3D &norm) {
  r.Reflect(norm);
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector (assignment)
/// \param v: vector
inline Vec4D &Vec4D::Reset(const Vec4D &v) { return *this = v; }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reset the vector components
/// \param _r0: time component
/// \param r1: space components in x-direction
/// \param r2: space components in y-direction
/// \param r3: space components in z-direction
inline Vec4D &Vec4D::Reset(VEC3D_T _r0, VEC3D_T r1, VEC3D_T r2, VEC3D_T r3) {
  r0 = _r0;
  r.x = r1;
  r.y = r2;
  r.z = r3;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator =
inline Vec4D &Vec4D::operator=(const Vec4D &v) {
  r0 = v.r0;
  r = v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator +=
inline Vec4D &Vec4D::operator+=(const Vec4D &v) {
  r0 += v.r0;
  r += v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator -=
inline Vec4D &Vec4D::operator-=(const Vec4D &v) {
  r0 -= v.r0;
  r -= v.r;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator *= scalar multiplication
/// \param d: scalar
inline Vec4D &Vec4D::operator*=(VEC3D_T d) {
  r0 *= d;
  r *= d;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Assignment operator *= scalar division
/// \param d: scalar
inline Vec4D &Vec4D::operator/=(VEC3D_T d) {
  r0 /= d;
  r /= d;
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator +
inline Vec4D Vec4D::operator+(const Vec4D &v) const {
  Vec4D temp(*this);
  return temp += v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator -
inline Vec4D Vec4D::operator-(const Vec4D &v) const {
  Vec4D temp(*this);
  return temp -= v;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator *
inline VEC3D_T Vec4D::operator*(const Vec4D &v) const {
  return r0 * v.r0 - r * v.r;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for left scalar multiplication
inline Vec4D operator*(VEC3D_T a, const Vec4D &v) {
  return Vec4D(a * v.r0, a * v.r);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar multiplication
inline Vec4D operator*(const Vec4D &v, VEC3D_T a) {
  return Vec4D(v.r0 * a, v.r * a);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator * for right scalar division
inline Vec4D operator/(const Vec4D &v, VEC3D_T a) {
  return Vec4D(v.r0 / a, v.r / a);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Negation
inline Vec4D operator-(const Vec4D &v) { return Vec4D(-v.r0, -v.r); }

}; // namespace IGCascade

namespace std {
inline ostream &operator<<(ostream &stream, const IGCascade::Vec4D &v);
inline istream &operator>>(istream &stream, IGCascade::Vec4D &v);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream insertion
inline ostream &operator<<(ostream &stream, const IGCascade::Vec4D &v) {
  v.DumpShort(stream);
  return stream;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Stream extraction
inline istream &operator>>(istream &stream, IGCascade::Vec4D &v) {
  char c;
  stream >> c >> v.r0 >> c >> v.r.x >> c >> v.r.y >> c >> v.r.z >> c;
  return stream;
}
} // namespace std

#endif // PHYSICS_VEC4D_H
