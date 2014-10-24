//-*-mode:c++; mode:font-lock;-*-

/*! \file vec3D.cpp
  vec3D class implementation file

  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n
  \author   Maciej Nicewicz           \n
            UCLA                      \n
	    nicewicz@physics.ucla.edu \n
  \author   Vladimir Vassiliev        \n
            UCLA                      \n
	    vvv@astro.ucla.edu        \n
  \author	Yusef Shafi
			UCLA
		yshafi@ucla.edu

  \date     7/25/2006
  \version  1.3
  \note
		1.3 changes all doubles to double-double precision using
		dd lib. Also, we add a uniform spherical direction function.
*/

/*! \example exmp_vec3D.cpp
    This is an example of how to use this class
 */

#include<limits>

#include "Vec3D.hpp"
//#include "Constants.hpp"
//#include "RandomNumbers.hpp"

using namespace IGCascade;

const double SMALLEST_ROTATION_ANGLE = 1.e-12;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Method for vector rotation around axis vector.
/// \param axis: axis of rotation vector [rad]
/*! \note
  The modulus of rotation angle is equal to the axis
  norm and is given in radians. The rotation of e_x
  around e_z by PI/2 is equal to e_y.
*/
void Vec3D::Rotate( const Vec3D& axis ) {

  VEC3D_T angle = axis.Norm();

  if(angle < SMALLEST_ROTATION_ANGLE)return;

  Vec3D e_axis = axis/angle;

  VEC3D_T par   = (*this)*e_axis;
  Vec3D p_axis = (*this)-par*e_axis;

#if 0
  VEC3D_T per   = p_axis.Norm();

  if(per == 0.0) return ;
  else p_axis=p_axis/per;

  *this=par*e_axis+cos(angle)*per*p_axis+sin(angle)*per*(e_axis^p_axis);
#else
  // "Rotation Formula" -- see Goldstein chap 4
  *this = par*e_axis + cos(angle)*p_axis + sin(angle)*(e_axis^p_axis);
#endif

  return;
}

/*void Vec3D::ScatterDirection(VEC3D_T dispersion, RandomNumbers& rng)
{
  VEC3D_T t = Norm();
  if(t==0)return;

  VEC3D_T phi = rng.Uniform() * Constants::num_pi;
  VEC3D_T theta = atan(dispersion*rng.Normal());

  Vec3D tangent_a;

  if((fabs(x)<=fabs(y))&&(fabs(x)<=fabs(z)))
	tangent_a = (*this)^Vec3D(1,0,0);
  else if(fabs(y)<=fabs(z))
	tangent_a = (*this)^Vec3D(0,1,0);
  else 
	tangent_a = (*this)^Vec3D(0,0,1);
  
  tangent_a /= tangent_a.Norm();
  
  Vec3D tangent_b((*this)^tangent_a);
  tangent_b /= tangent_b.Norm();

#if 1
  *this = (*this)*cos(theta) + (tangent_a*cos(phi)+
  tangent_b*sin(phi))*t*sin(theta);
#else  
  Vec3D axis = tangent_a*cos(phi) + tangent_b*sin(phi);
  Rotate(axis*theta);
#endif
}*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator &= of rotation composition
/*! \note 
  Composition is opposite to the sense of matrix
  multiplication. The composition r=r1&r2 is equivalent to 
  a rotation first by r1 then by r2.

  For example, for two rotations r1 and r2, and for any vector p

  Vec3D r1(?,?,?);
  Vec3D r2(?,?,?);
  Vec3D p(?,?,?);

  Vec3D r=r1&r2;

  p.Rotate(r1);   \/  is equivalent to p.Rotate(r)
  p.Rotate(r2);   /\ 
*/

/* 
   Composition algorithm comes from consideration of rotation
   as 2-spinors... i.e. write the both rotations in terms of Pauli
   matrices, multiply and then compare terms. See any QM book or
   Goldstein chap 4.

   $\mathbold{R_1} = \mathbold{1} \cos\theta_1 / 2 + 
    \hat{n}_1 . \vec{\mathbold{sigma}} \sin\theta_1 / 2 $

   $\mathbold{R_2} = \mathbold{1} \cos\theta_2 / 2 + 
    \hat{n}_2 . \vec{\mathbold{sigma}} \sin\theta_2 / 2 $

    Multiply the matrices using the following, collect terms
    and compare with $R_1$ or $R_2$ to deduce the composite
    angle and axis of rotation.

   $(\vec{\mathbold{sigma}}.\vec{A})
    (\vec{\mathbold{sigma}}.\vec{B}) = \mathbold{1}\vec{A}.\vec{B}+
    i\vec{\mathbold{sigma}}(\vec{A}\cross\vec{B})$
*/

Vec3D& Vec3D::operator &= (const Vec3D& r)
{
#define R1 (*this)
#define R2 r

  VEC3D_T r1_theta = R1.Norm();
  VEC3D_T r2_theta = R2.Norm();

  if(r1_theta == 0.0)
    {
      // R1 is zero so set to R2
      *this = R2;
    }
  else if(r2_theta == 0.0)
    {
      // R2 is zero so set to R1
      *this = R1;
    }
  else
    {
      VEC3D_T sin_half_r1 = sin(r1_theta/2.0);
      VEC3D_T cos_half_r1 = cos(r1_theta/2.0);
      VEC3D_T sin_half_r2 = sin(r2_theta/2.0);
      VEC3D_T cos_half_r2 = cos(r2_theta/2.0);

      Vec3D r1_hat = R1/r1_theta;
      Vec3D r2_hat = R2/r2_theta;

      VEC3D_T coshalftheta =
	cos_half_r1*cos_half_r2-
	sin_half_r1*sin_half_r2*(r1_hat*r2_hat);

      Vec3D axis(r1_hat*sin_half_r1*cos_half_r2+
		 r2_hat*sin_half_r2*cos_half_r1-
		 (r1_hat^r2_hat)*sin_half_r1*sin_half_r2);

      VEC3D_T sinhalftheta = axis.Norm();

      VEC3D_T halftheta=atan(sinhalftheta/coshalftheta);

      *this = axis/sinhalftheta*halftheta*2.0;
    }
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components and vector norm in square
void Vec3D::Dump( std::ostream& stream ) const
{
  stream << std::endl;
  stream << " .X:   " << x << std::endl;
  stream << " .Y:   " << y << std::endl;
  stream << " .Z:   " << z << std::endl;
  stream << "Norm2: " << Norm2() << std::endl;
}
/*ORIGINAL
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components in short form
void Vec3D::DumpShort( std::ostream& stream ) const
{
  stream << "( " << x << ' ' << y << ' ' << z << " )";
}
*/

//EDIT 6/29/06
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Print vector components in short form
void Vec3D::DumpShort( std::ostream& stream ) const
{
  stream << x << ',' << y << ',' << z;
}

/*
#ifdef TESTMAIN

// Compile with: g++ -O3 -DTESTMAIN -o test Vec3D.cpp -lm

#include<sstream>
#include<vector>

#include"RandomNumbers.hpp"

int main()
{
  RandomNumbers rng("random.seeds");

  Vec3D a(10,20,30);
  Vec3D b(1,2,3);
  std::cout << "sizeof(Vec3D): " << sizeof(Vec3D) << std::endl;
  std::cout << "sizeof(Vec3D*): " << sizeof(Vec3D*) << std::endl << std::endl;

  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Norm2():" << a.Norm2() << "  "
	    << "b.Norm2():" << b.Norm2() << std::endl;
  std::cout << "a.Norm():" << a.Norm() << "  "
	    << "b.Norm():" << b.Norm() << std::endl;

  std::cout << "a+=b: " << (a+=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a-=b: " << (a-=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;

  std::cout << "a*=1.5: " << (a*=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a/=1.5: " << (a/=1.5) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(): " << a.Reset() << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Reset(b): " << a.Reset(b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(0,1,0): " << a.Reset(0,1,0) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(0,0,1): " << a.Reset(0,0,1) << std::endl;
  std::cout << "b.Reset(0,0,pi/4): " << b.Reset(0,0,M_PI/2) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl;
  std::cout << "a.Rotate(b)" << std::endl;
  a.Rotate(b);
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,-2,3): " << a.Reset(1,-2,3) << std::endl;
  std::cout << "a.P()" << std::endl;
  a.P();
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(1,0,0): " << b.Reset(1,0,0) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,1,0): " << b.Reset(0,1,0) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  std::cout << "a.Reset(1,0,0): " << a.Reset(1,0,0) << std::endl;
  std::cout << "b.Reset(0,0,1): " << b.Reset(0,0,1) << std::endl;
  std::cout << "a^=b:" << (a^=b) << std::endl;
  std::cout << "a:" << a << " b:" << b << std::endl << std::endl;

  Vec3D c;

  std::cout << "c=b=a=():" << (c=b=a=Vec3D()) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "c=b=(a=())+(1,0,0)):"
	    << (c=b=(a=Vec3D())+Vec3D(1,0,0)) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "a=(1,0,0):" << (a=Vec3D(1,0,0)) << std::endl;
  std::cout << "b=(1,-1,0):" << (b=Vec3D(1,-1,0)) << std::endl;
  std::cout << "c=a+b:" << (c=a+b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "c=a-b:" << (c=a-b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "a*b:" << (a*b) << std::endl;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  std::cout << "(3,4,5)*(2,0,0): " << (Vec3D(3,4,5)*Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)*(0,3,0): " << (Vec3D(3,4,5)*Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)*(0,0,4): " << (Vec3D(3,4,5)*Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)*(2,3,0): " << (Vec3D(3,4,5)*Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)*(2,0,4): " << (Vec3D(3,4,5)*Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)*(0,3,4): " << (Vec3D(3,4,5)*Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)*(2,3,4): " << (Vec3D(3,4,5)*Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)^(2,0,0): " << (Vec3D(3,4,5)^Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)^(0,3,0): " << (Vec3D(3,4,5)^Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)^(0,0,4): " << (Vec3D(3,4,5)^Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)^(2,3,0): " << (Vec3D(3,4,5)^Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)^(2,0,4): " << (Vec3D(3,4,5)^Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)^(0,3,4): " << (Vec3D(3,4,5)^Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)^(2,3,4): " << (Vec3D(3,4,5)^Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)+(2,0,0): " << (Vec3D(3,4,5)+Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)+(0,3,0): " << (Vec3D(3,4,5)+Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)+(0,0,4): " << (Vec3D(3,4,5)+Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)+(2,3,0): " << (Vec3D(3,4,5)+Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)+(2,0,4): " << (Vec3D(3,4,5)+Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)+(0,3,4): " << (Vec3D(3,4,5)+Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)+(2,3,4): " << (Vec3D(3,4,5)+Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "(3,4,5)-(2,0,0): " << (Vec3D(3,4,5)-Vec3D(2,0,0)) << std::endl;
  std::cout << "(3,4,5)-(0,3,0): " << (Vec3D(3,4,5)-Vec3D(0,3,0)) << std::endl;
  std::cout << "(3,4,5)-(0,0,4): " << (Vec3D(3,4,5)-Vec3D(0,0,4)) << std::endl;
  std::cout << "(3,4,5)-(2,3,0): " << (Vec3D(3,4,5)-Vec3D(2,3,0)) << std::endl;
  std::cout << "(3,4,5)-(2,0,4): " << (Vec3D(3,4,5)-Vec3D(2,0,4)) << std::endl;
  std::cout << "(3,4,5)-(0,3,4): " << (Vec3D(3,4,5)-Vec3D(0,3,4)) << std::endl;
  std::cout << "(3,4,5)-(2,3,4): " << (Vec3D(3,4,5)-Vec3D(2,3,4)) << std::endl;
  std::cout << std::endl;

  std::cout << "-(3,4,5): " << (-Vec3D(3,4,5)) << std::endl;
  std::cout << "1.5*(3,4,5): " << (1.5*Vec3D(3,4,5)) << std::endl;
  std::cout << "(3,4,5)*1.5: " << (Vec3D(3,4,5)*1.5) << std::endl;
  std::cout << "(3,4,5)/1.5: " << (Vec3D(3,4,5)/1.5) << std::endl << std::endl;

  std::cout << "std::istringstream s(\"(10,11,12)\"); s >> a;" << std::endl;
  std::istringstream s("(10,11,12)");
  s >> a;
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  for(unsigned nrot = 1; nrot<=20; nrot++)
    {
      VEC3D_T mean_residual = 0;
      VEC3D_T max_residual = 0;

      for(unsigned n=0; n<1000; n++)
	{
	  VEC3D_T phi = rng.Uniform() * M_PI * 2.0;
	  VEC3D_T costheta = rng.Uniform() * 2.0 - 1.0;
	  VEC3D_T sintheta = sin(acos(costheta));
	  Vec3D a(sintheta*cos(phi),sintheta*sin(phi),costheta);
	  Vec3D b(a);
	  Vec3D rc;

	  std::vector<Vec3D> R;
	  std::vector<Vec3D> Rc;

	  for(unsigned i=0; i<nrot; i++)
	    {
	      phi = rng.Uniform() * M_PI * 2.0;
	      costheta = rng.Uniform() * 2.0 - 1.0;
	      sintheta = sin(acos(costheta));
	      Vec3D r(sintheta*cos(phi),sintheta*sin(phi),costheta);
	      r *= rng.Uniform()*2.0*M_PI;

	      b.Rotate(r);
	      rc &= r;

	      R.push_back(r);
	      Rc.push_back(rc);
	    }

	  a.Rotate(rc);

	  Vec3D c(a-b);
	  VEC3D_T residual = c.Norm();
	  if(residual>max_residual)max_residual=residual;
	  mean_residual+=residual;

#if 0
	  if(residual>1e-14)
	    {
	      for(unsigned i=0; i<R.size(); i++)
		std::cout << R[i].Norm()
			  << ' '  << R[i] << ' ' << Rc[i] << std::endl;
	      std::cout << a << ' ' << b << ' ' << c << std::endl;
	    }
#endif
	}

      std::cout << nrot << " rotations: mean residual: "
		<< mean_residual/1000.0
		<< " max residual: " << max_residual << std::endl;
    }

  Vec3D r1=M_PI/2*Vec3D(1,0,0);
  Vec3D r2=M_PI/2*Vec3D(0,1,0);
  Vec3D rc = r1 & r2;
  std::cout << "r1:" << r1 << " r2:" << r2 << " rc:" << rc << std::endl;

  a=Vec3D(1,0,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(0,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(0,0,2); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;

  a=Vec3D(1,1,0); b=a; b.Rotate(r1); b.Rotate(r2); c=a; c.Rotate(rc);
  std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl << std::endl;

  a=Vec3D(.1,0,0);
  a &= a & a & a;
  std::cout << "a:" << a << std::endl;
}

#endif // TESTMAIN
*/
