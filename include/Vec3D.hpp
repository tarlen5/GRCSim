//-*-mode:c++; mode:font-lock;-*-

/*! \file Vec3D.hpp
  Vec3D class header file
  
  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n
  \author   Maciej Nicewicz           \n
            UCLA                      \n
	    nicewicz@physics.ucla.edu \n
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n
  \author   Yusef Shafi
            UCLA
	    yshafi@ucla.edu

  \date     7/25/2006
  \version  1.3
  \note
	    1.3 changes all doubles to double-double precision using
	        dd lib. Also, we add a uniform spherical direction function.
*/

#ifndef IGCASCADE_VEC3D_H
#define IGCASCADE_VEC3D_H

#ifndef __VS_NO_IOSTREAM
#include <iostream>
#endif

#include<qd/dd_real.h>
//#include<qd/dd.h>
//#include<qd/qd.h>

//#include "RandomNumbers.hpp"
#include "TRandom3.h"

#include <cmath>

#define VEC3D_PI dd_real::_pi

typedef dd_real VEC3D_T;

/*!  \class Vec3D
     \brief 3 dimensional vector class

     This class defines 3D vectors, set of operations 
     in vector field, scalar and vector products, 
     as well as rotation, parity transformation, negation,
     and normalization.

*/

//class RandomNumbers;

namespace IGCascade
{
  class Vec3D
  {
  public:
    inline Vec3D();                                 //!<default constructor
    inline Vec3D(const Vec3D& v);                   //!<copy constructor
    inline Vec3D(VEC3D_T _x,VEC3D_T _y,VEC3D_T _z); //!<overloaded constructor

    inline VEC3D_T Norm() const;     //!<calculates norm
    inline VEC3D_T Norm2() const;    //!<calculates scalar product with itself
    void Rotate(const Vec3D& axis);  //!<rotates vector around axis
    inline void P();                 //!<parity transformation
    inline void Reflect(const Vec3D& norm);  //!<reflect in normal
    static inline Vec3D UniformSphereDirection(TRandom3* rng); 
                           //!<choose a (normalized) direction randomly

    //void ScatterDirection(VEC3D_T dispersion, RandomNumbers& rng);

    inline Vec3D& Reset(const Vec3D& v = Vec3D());
    inline Vec3D& Reset(VEC3D_T _x, VEC3D_T _y, VEC3D_T _z);

    inline Vec3D& operator = ( const Vec3D& v );  //!<assignment
    inline Vec3D& operator += ( const Vec3D& v ); //!<assignment: addition
    inline Vec3D& operator -= ( const Vec3D& v ); //!<assignment: subtraction 
    inline Vec3D& operator ^= ( const Vec3D& v ); //!<vector product  
    inline Vec3D& operator *= (VEC3D_T d); //!<assignment: multiply by scaler
    inline Vec3D& operator /= (VEC3D_T d); //!<assignment: divide by scaler

    Vec3D& operator &=(const Vec3D& r);//!<assignment: composition of rotations

    inline Vec3D  operator + (const Vec3D& v) const;  //!<addition
    inline Vec3D  operator - (const Vec3D& v) const;  //!<subtraction
    inline VEC3D_T operator * (const Vec3D& v) const;  //!<scalar product
    inline Vec3D  operator ^ (const Vec3D& v) const;  //!<vector product

    inline Vec3D  operator & (const Vec3D& v) const;  //!<addition of rotations

#ifndef __VS_NO_IOSTREAM
    void Dump(std::ostream& stream = std::cout) const; //!<prints coordinates
    void DumpShort(std::ostream& stream = std::cout) const; //!<prints coords
#endif

    VEC3D_T x, y, z;   //!<components

  private:

  };

  inline Vec3D operator - ( const Vec3D& v );           //!<negation
  inline Vec3D operator * ( VEC3D_T d, const Vec3D& v ); //!<left scalar mult.
  inline Vec3D operator * ( const Vec3D& v, VEC3D_T d ); //!<right scalar mult.
  inline Vec3D operator / (const Vec3D& v,VEC3D_T d); //!<right scalar division

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Default class constructor
  inline Vec3D::Vec3D(): x(), y(), z()
  {
    // nothing to see here
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Copy constructor
  // \param o: Vector to copy
  inline Vec3D::Vec3D(const Vec3D& o): x(o.x), y(o.y), z(o.z)
  {
    // nothing to see here
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Overloaded class constructor
  /// \param _x: x
  /// \param _y: y
  /// \param _z: z
  inline Vec3D::Vec3D( VEC3D_T _x, VEC3D_T _y, VEC3D_T _z ): 
    x(_x), y(_y), z(_z)
  {
    // nothing to see here
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Method for vector norm
  inline VEC3D_T Vec3D::Norm() const
  {
    return sqrt(x*x + y*y + z*z);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Method for vector norm in square
  inline VEC3D_T Vec3D::Norm2() const
  {
    return x*x + y*y + z*z;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Method for parity transformation of vector
  inline void Vec3D::P()
  {
    *this=-(*this);
    return;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Method for reflection of vector in mirror
  inline void Vec3D::Reflect(const Vec3D& norm)  //!reflect in normal
  {
    VEC3D_T n2 = norm.Norm2();
    if(n2 == (VEC3D_T)(0.0))return;
    *this -= norm * ((VEC3D_T)(2.0)*((*this)*norm)/n2);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Reset the vector (assignment)
  /// \param v: vector
  inline Vec3D& Vec3D::Reset(const Vec3D& v)
  {
    return *this = v;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Reset the vector components
  /// \param _x: x
  /// \param _y: y
  /// \param _z: z
  inline Vec3D& Vec3D::Reset(VEC3D_T _x, VEC3D_T _y, VEC3D_T _z)
  {
    x = _x;
    y = _y;
    z = _z;
    return *this;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator =
  inline Vec3D& Vec3D::operator = ( const Vec3D& v )
  {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator +=
  inline Vec3D& Vec3D::operator += ( const Vec3D& v )
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator -=
  inline Vec3D& Vec3D::operator -= ( const Vec3D& v )
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator ^= vector multiplication
  /// \param v: vector
  /*! \note 
    Normally ^ operator is a bitwise exclusive OR, which 
    has precedence lower than * / and even lower than + -.
    Thus it is executed last if no brackets are used.
    Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
    See examples file exmp_Vec3D.cpp for additional
    comment(s)
  */
  inline Vec3D& Vec3D::operator ^= ( const Vec3D& v )
  {
    Vec3D temp(y*v.z - z*v.y,
	       z*v.x - x*v.z,
	       x*v.y - y*v.x);
    *this=temp;
    return *this;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator *= scalar multiplication
  /// \param d: scalar
  inline Vec3D& Vec3D::operator *= ( VEC3D_T d )
  {
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Assignment operator /= scalar division
  /// \param d: scalar
  inline Vec3D& Vec3D::operator /= ( VEC3D_T d )
  {
    x /= d;
    y /= d;
    z /= d;
    return *this;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator +
  inline Vec3D Vec3D::operator + ( const Vec3D& v ) const
  {
    Vec3D temp(*this);
    return temp += v;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator -
  inline Vec3D Vec3D::operator - ( const Vec3D& v ) const
  {
    Vec3D temp(*this);
    return temp -= v;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator *
  inline VEC3D_T Vec3D::operator * ( const Vec3D& v ) const
  {
    return x*v.x + y*v.y + z*v.z;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator ^ of vector product
  /*! \note 
    Normally ^ operator is a bitwise exclusive OR, which 
    has precedence lower than * / and even lower than + -.
    Thus it is executed last if no brackets are used.
    Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
    See examples file exmp_Vec3D.cpp for additional
    comment(s)
  */
  inline Vec3D Vec3D::operator ^ ( const Vec3D& v ) const
  {
    Vec3D temp(*this);
    return temp ^= v;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator & of rotation composition
  /*! \note 
    Composition is opposite to the sense of matrix
    multiplication. If r=r1&r2 then rotation r1 happens
    first, followed by r2

    Normally & operator is a bitwise exclusive AND, which 
    has precedence lower than * / and even lower than + -.
    Thus it is executed last if no brackets are used.
    Exmp: r=7.5*a+l*b-c*m+2*b^c=(7.5*a+l*b-c*m+2*b)^c
    See examples file exmp_Vec3D.cpp for additional
    comment(s)
  */
  inline Vec3D Vec3D::operator & ( const Vec3D& v ) const
  {
    Vec3D temp(*this);
    return temp &= v;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator * for left scalar multiplication
  inline Vec3D operator * ( VEC3D_T d, const Vec3D& v )
  {
    Vec3D temp(d*v.x, d*v.y, d*v.z );
    return temp;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator * for right scalar multiplication
  inline Vec3D operator * ( const Vec3D& v, VEC3D_T d )
  {
    Vec3D temp( v.x*d, v.y*d, v.z*d );
    return temp;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Operator * for right scalar division
  inline Vec3D operator / ( const Vec3D& v, VEC3D_T d )
  {
    Vec3D temp( v.x/d, v.y/d, v.z/d );
    return temp;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Negation
  inline Vec3D operator - ( const Vec3D& v )
  {
    Vec3D temp( -v.x, -v.y, -v.z );
    return temp;
  }
  
  inline Vec3D Vec3D::UniformSphereDirection(TRandom3* rng)
  {
    VEC3D_T phi = (2.*VEC3D_PI)*rng->Uniform();
    VEC3D_T cos_theta = 1. - 2.*rng->Uniform();
    
    VEC3D_T sin_theta = sqrt(1.-cos_theta*cos_theta);
    
    VEC3D_T x,y,z;
    z = cos_theta;
    x = sin_theta*cos(phi);
    y = sin_theta*sin(phi);
    
    Vec3D v(x,y,z);
    
    return v;
  }

} // namespace IGCascade

#ifndef __VS_NO_IOSTREAM

namespace std
{
  inline ostream& operator << (ostream& stream, const IGCascade::Vec3D& v);
  inline istream& operator >> (istream& stream, IGCascade::Vec3D& v);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Stream insertion
  inline ostream& operator << (ostream& stream, const IGCascade::Vec3D& v)
  {
    v.DumpShort(stream);
    return stream;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Stream extraction
  inline istream& operator >> (istream& stream, IGCascade::Vec3D& v)
  {
    char c;
    stream >> c >> v.x >> v.y >> v.z >> c;
    return stream;
  }
}

#endif

#endif // IGCascade_VEC3D_H
