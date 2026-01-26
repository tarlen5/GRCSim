/*!
  \file    IGMFPropagator.hpp
           Abstract base class for intergalactic magnetic field propagators.

  \author  Timothy C. Arlen

  \date    01/26/2026

  \note    Provides common interface for magnetic field propagation models.
           Generated for refactoring based on MagneticGrid and
           MFTurbulentContinuous.

*/

#ifndef IGCASCADE_IGMFPROPAGATOR_H
#define IGCASCADE_IGMFPROPAGATOR_H

#include "HighPrecProp.hpp"
#include "PhysicsConstants.hpp"
#include "RandomNumbers.hpp"
#include "RelParticle.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "convert.hpp"

namespace IGCascade {

class IGMFPropagator {
protected:
  RandomNumbers *m_rng;

public:
  inline IGMFPropagator(RandomNumbers *rng, VEC3D_T bmag)
      : m_rng(rng), m_bmag(bmag), m_DE("1.0E-25") {}
  virtual ~IGMFPropagator() = default;

  virtual void PropagateBFieldRedshift(
      RelParticle &photon, RelParticle *&lepton, Vec3D &n_eo,
      VEC3D_T &prop_length, VEC3D_T &delta_z
  ) = 0;

  VEC3D_T m_bmag;
  VEC3D_T m_DE;
};

} // namespace IGCascade

#endif // IGCASCADE_IGMFPROPAGATOR_H