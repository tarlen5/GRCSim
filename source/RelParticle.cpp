/*! \file RelParticle.cpp
  RelParticle class implementation file

  \author   Yusef Shafi             \n
            UCLA                      \n
      yshafi@ucla.edu     \n

  \author   Timothy C. Arlen    \n
            Department of Physics and Astronomy       \n
            UCLA                                      \n
      email:   arlen@astro.ucla.edu             \n

  \date     08/07/2006

  \version  1.1
            1.2 - 08/15/08

  \revision 08/15/08 - Updated and fixed Dump() function

  \note
*/

#include "RelParticle.hpp"

using namespace IGCascade;

// const double m_pi = 3.14159265358979;

namespace IGCascade {
Vec4D RelParticle::GetUcm(std::vector<RelParticle> &particles) {
  Vec3D momenta;
  VEC3D_T mass = 0;
  for (unsigned ctr = 0; ctr < particles.size(); ctr++) {
    momenta += particles[ctr].m_p4.r;
    mass += particles[ctr].m_p4.r0;
  }

  Vec3D u_cm;

  u_cm = momenta / mass;

  VEC3D_T gamma = 1. / sqrt(1. - (u_cm * u_cm));

  Vec4D u_cm4(gamma, gamma * u_cm);

  return u_cm4;
}

void RelParticle::BoostSystem(Vec4D &u4, std::vector<RelParticle> &particles) {
  for (unsigned ctr = 0; ctr < particles.size(); ctr++) {
    particles[ctr].Boost(u4);
  }
}

/*Vec3D RelParticle::UniformSphereDirection()
  {
  RandomNumbers rng("uniform_dist.seeds");
  VEC3D_T phi = (2.*m_pi)*rng.Uniform();
  VEC3D_T u = -1. + 2.*rng.Uniform();

  VEC3D_T x,y,z;
  z = u;
  x = sqrt(1.-u*u)*cos(phi);
  y = sqrt(1.-u*u)*sin(phi);

  Vec3D v(x,y,z);

  return v;
  }*/

} // namespace IGCascade

void RelParticle::Dump(std::ostream &stream) const {
  stream << std::endl;
  stream << "*****Particle*****" << std::endl;
  // stream << " Type: " << tag << std::endl;
  stream << " Charge:	" << m_q << std::endl;
  stream << " Tag:	" << m_tag << std::endl;
  stream << " Rest mass: " << m_m0 << std::endl;
  stream << " redshift: " << m_z << std::endl;
  stream << " redshift_int: " << m_z_int << std::endl;
  stream << " redshift_s: " << m_z_s << std::endl;
  stream << " redshift_s_int: " << m_z_s_int << std::endl;

  stream << std::endl;
  stream << " Position:" << std::endl;
  stream << " .T:   " << m_r4.r0 << std::endl;
  stream << " .X:   " << m_r4.r.x << std::endl;
  stream << " .Y:   " << m_r4.r.y << std::endl;
  stream << " .Z:   " << m_r4.r.z << std::endl;
  // stream << " Norm2: " << m_r4.Norm2() << std::endl;

  stream << std::endl;
  stream << " Energy-Momentum:" << std::endl;
  stream << " .Energy (E):   " << m_p4.r0 << std::endl;
  stream << " .P_x:   " << m_p4.r.x << std::endl;
  stream << " .P_y:   " << m_p4.r.y << std::endl;
  stream << " .P_z:   " << m_p4.r.z << std::endl;
  // stream << " Norm2: " << m_p4.Norm2() << std::endl;

  // stream << std::endl;
  stream << "****/Particle*****" << std::endl;
  stream << std::endl;
}

void RelParticle::DumpShort(std::ostream &stream) const {
  stream << m_r4.r0 << ", " << m_r4.r.x << ", " << m_r4.r.y << ", " << m_r4.r.z;
}
