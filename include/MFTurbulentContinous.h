/*!
  \file    MFTurbulentContinous.h
           Header file for MFTurbulentContinous.cpp

  \author  Tim Arlen
           timothyarlen@gmail.com

  \date    01/17/2026

  \note    This class implements a continuous turbulent magnetic field
           model for intergalactic magnetic fields, removes need for a
           grid-based approach, any file locking, I/O, etc. by using a
           procedural noise functions like Perlin which is designed to
           generate such a pseudo-random, continuous field.

*/

#ifndef IGCASCADE_MFTURBULENTCONTINOUS_H
#define IGCASCADE_MFTURBULENTCONTINOUS_H

#include <cstdint>
#include <random>

#include "FastNoiseLite.h"
#include "HighPrecProp.hpp"
#include "PhysicsConstants.hpp"
#include "RandomNumbers.hpp"
#include "RelParticle.hpp"
#include "Vec3D.hpp"
#include "convert.hpp"

namespace IGCascade {

class MFTurbulentContinous {
public:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors///////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  MFTurbulentContinous(RandomNumbers *rng, VEC3D_T b_mag, VEC3D_T coh_len);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Member Functions///////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void PropagateBFieldRedshift(
      RelParticle &photon, RelParticle *&lepton, Vec3D &n_eo,
      VEC3D_T &prop_length, VEC3D_T &delta_z
  );

  Vec3D GetMagneticFieldDirectionAtPosition(const Vec3D &position);

  inline VEC3D_T getCoherenceLength() { return m_coh_len; }
  inline VEC3D_T getBFieldMagnitude() { return m_bmag; }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Public Data Members////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  VEC3D_T m_DE; // minimum precision for root finding at z_s = 0 surface

private:
  Vec3D GetUpdatedPositionWithMFDeflection(
      RelParticle *&lepton, VEC3D_T &prop_length_step, const Vec3D &n_eo,
      const Vec3D &b_field_dir
  );

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Private Data Members////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  RandomNumbers *m_rng;
  VEC3D_T m_bmag;    // B_magnitude [gauss]
  VEC3D_T m_coh_len; // Coherence Length [cm]
  uint32_t m_seed;   // Seed for noise generation

  // Noise generators for each component
  FastNoiseLite m_noise_x, m_noise_y, m_noise_z;
};
} // namespace IGCascade

#endif // IGCASCADE_MFTURBULENTCONTINOUS_H