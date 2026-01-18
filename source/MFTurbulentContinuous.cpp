/*!
-------------------------------------------------------------------------------
  \file     MFTurbulentContinous.cpp
            Propagates charged leptons and photons through magnetic fields
            once the propagation length is defined.

  \author   Timothy C. Arlen
            timothyarlen@gmail.com

  \date     01/18/2026
-------------------------------------------------------------------------------
*/

#include "MFTurbulentContinuous.h"

namespace IGCascade {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MFTurbulentContinuous::MFTurbulentContinuous(
    RandomNumbers *rng, VEC3D_T b_mag, VEC3D_T coh_len
) {

  std::cout << "Initializing MFTurbulentContinuous with B = " << Double(b_mag)
            << " G, Coherence Length = " << Double(coh_len) << " Mpc."
            << std::endl;

  m_rng = rng;
  m_bmag = b_mag;      // [gauss]
  m_coh_len = coh_len; // [Mpc]
  m_DE = "1.0E-25";

  if (rng->IsSeedSet()) {
    m_seed = rng->GetSeed();
  } else {
    std::random_device rd;
    m_seed = rd();
  }

  m_noise_x.SetSeed(m_seed);
  m_noise_y.SetSeed(m_seed + 1); // Offset for decorrelation
  m_noise_z.SetSeed(m_seed + 2);

  double d_coh_len = Double(m_coh_len);
  float freq = 1.0f / static_cast<float>(d_coh_len); // Units: e.g., Mpc^{-1}
  m_noise_x.SetFrequency(freq);
  m_noise_y.SetFrequency(freq);
  m_noise_z.SetFrequency(freq);
  m_noise_x.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
  m_noise_y.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
  m_noise_z.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
}

void MFTurbulentContinuous::PropagateBFieldRedshift(
    RelParticle &photon, RelParticle *&lepton, Vec3D &n_eo,
    VEC3D_T &prop_length, VEC3D_T &delta_z
)
/*
Once a propagation length (and corresponding delta_z depending on starting
redshift) has been determined for the lepton (and it's associated scattered
photon), this function propagates both particles through the magnetic field,
translating and rotating them accordingly and directly modifying their
4-positions and 4-momenta in the RelParticle class.

The most involved aspect of this class is making sure the propagation doesn't
take you past the z=0 surface, in which case the propagatoin length is
shortened to end exactly at z=0. Secondly, the propagation through constant
magnetic field occurs until the propagation length exceeds the coherence
length, at which point, a new magnetic field vector is queried at the new
position and the rest of the propagation length is done in that new field, etc.

  \param photon      - The photon which will be rotated and translated

  \param lepton      - The lepton which will be rotated and translated

  \param n_eo        - Lepton direction before IC scat with photon.

  \param prop_length - Propagation length of lepton                     [cm]

  \param delta_z     - change in redshift of the lepton over the
                       distance of the propagation length.
*/
{

  // Safety checks
  if (!lepton) {
    std::cerr << "Error: lepton pointer is null!" << std::endl;
    std::abort();
  }
  // Example checks for lepton members (customize as needed)
  // If m_r4 or m_p4 are pointers, check them too:
  // if (!lepton->m_r4) { ... }
  // if (!lepton->m_p4) { ... }
  /*
  break source/MFTurbulentContinuous.cpp:134
  set args 3000 1e-17 1 0.14 1 10
  */

  // Check for zero-length direction vector before normalization
  Vec3D test_b_field = GetMagneticFieldDirectionAtPosition(lepton->m_r4.r);
  if (test_b_field.Norm() == 0.0) {
    std::cerr
        << "Error: Magnetic field direction vector has zero length at position "
        << lepton->m_r4.r << std::endl;
    std::abort();
  }

  // cout<<"  Starting lepton propagation through turbulent MF..."<<endl;

  VEC3D_T prop_length_remaining = prop_length;
  VEC3D_T coherence_length_cm = PhysConst::MPC_TO_CM * m_coh_len;

  VEC3D_T distance_rtol = "1.0E-12";
  VEC3D_T prop_length_tol = distance_rtol * prop_length;

  while (prop_length_remaining > prop_length_tol) {
    VEC3D_T prop_length_next_step = prop_length_remaining;

    // check if probp length remaining is larger than coherence length
    // if it is, adjust prop length next step to coh length.
    // Then check whether prop_length_next_step would cross z=0 surface.
    // If it does, shorten prop length next step to end exactly at z=0 surface

    if (prop_length_remaining > coherence_length_cm) {
      prop_length_next_step = coherence_length_cm;
    }

    // Next calculate would-be next position after this prop_length_next_step
    // including deflection by magnetic field, and then check if it would
    // cross z=0 surface.
    Vec3D b_field_dir = GetMagneticFieldDirectionAtPosition(lepton->m_r4.r);
    // std::cout << "    b_field_dir at position " << lepton->m_r4.r
    //       << " = (" << b_field_dir.x << ", "
    //       << b_field_dir.y << ", "
    //       << b_field_dir.z << "), Norm = "
    //       << b_field_dir.Norm() << std::endl;
    Vec3D r_new_if_propagated = GetUpdatedPositionWithMFDeflection(
        lepton, prop_length_next_step, n_eo, b_field_dir
    );

    // Update prop_length_remaining here.
    VEC3D_T delta_z = GetLeptonDelta_z(prop_length_next_step, lepton);
    VEC3D_T delta_time = 0.0;
    VEC3D_T delta_zs = GetPhotonDelta_zs(
        prop_length_next_step, lepton, r_new_if_propagated, delta_time
    );

    // If z_s = 0 surface NOT crossed, do normal propagation and updating:
    if ((lepton->m_z_s - delta_zs) >= "0.0") {

      // Here is where the seg fault occurs. Why?
      PropagateConstantMF(
          photon, lepton, m_bmag, prop_length_next_step, r_new_if_propagated,
          delta_time, n_eo, delta_z, delta_zs, b_field_dir
      );
      prop_length_remaining -= prop_length_next_step;

    } else {
      // Now we will have crossed the z_s = 0 surface.
      // Need to shorten prop length next step to end exactly at z_s = 0

      VEC3D_T pl_right = prop_length_next_step;
      VEC3D_T pl_left = "0.0";
      Vec3D r_e = lepton->m_r4.r;

      // Bisection method to find the correct propagation length z_s=0
      while (fabs(lepton->m_z_s - delta_zs) > m_DE) {
        prop_length_next_step = (pl_right + pl_left) / 2.0;
        r_new_if_propagated = GetUpdatedPositionWithMFDeflection(
            lepton, prop_length_next_step, n_eo, b_field_dir
        );
        delta_zs = GetPhotonDelta_zs(
            prop_length_next_step, lepton, r_new_if_propagated, delta_time
        );

        if ((lepton->m_z_s - delta_zs) > 0.0) {
          pl_left = prop_length_next_step;
        } else {
          pl_right = prop_length_next_step;
        }
      }
      ///////////////////////////////////////////////////////////

      // Now calculate all values once prop_length_next_step to z_s=0 surface
      // has been computed:
      delta_z = GetLeptonDelta_z(prop_length_next_step, lepton);

      PropagateConstantMF(
          photon, lepton, m_bmag, prop_length_next_step, r_new_if_propagated,
          delta_time, n_eo, delta_z, delta_zs, b_field_dir
      );

      prop_length_remaining = "0.0"; // End propagation
    }

  } // end while prop length remaining > 0

  // cout<<"  Finished lepton propagation through turbulent MF."<<endl;
}

Vec3D MFTurbulentContinuous::GetUpdatedPositionWithMFDeflection(
    RelParticle *&lepton, VEC3D_T &prop_length_step, const Vec3D &n_eo,
    const Vec3D &b_field_dir
) {
  VEC3D_T ze = lepton->m_z;
  VEC3D_T p_o = lepton->m_p4.r.Norm() / (1.0 + ze);

  VEC3D_T Kappa_o = CGS_C * 1.E-8;
  Kappa_o *= lepton->m_q * m_bmag / p_o;

  Vec3D DeltaR = n_eo * prop_length_step;
  // The following gives the correction for magnetic deflection
  if (Kappa_o > "0.0") {

    VEC3D_T KappaL = Kappa_o * prop_length_step;
    Vec3D B = (n_eo - (n_eo * b_field_dir) * b_field_dir) *
              (sin(KappaL) - KappaL) / Kappa_o;
    Vec3D C = (n_eo ^ b_field_dir) * (1.0 - cos(KappaL)) / Kappa_o;
    DeltaR += (B + C);
  }

  Vec3D r_new = lepton->m_r4.r + DeltaR;

  return r_new;
}

Vec3D MFTurbulentContinuous::GetMagneticFieldDirectionAtPosition(
    const Vec3D &position
) {
  float x = static_cast<float>(Double(position.x));
  float y = static_cast<float>(Double(position.y));
  float z = static_cast<float>(Double(position.z));

  // Convert position to Mpc units:
  x /= static_cast<float>(Double(PhysConst::MPC_TO_CM));
  y /= static_cast<float>(Double(PhysConst::MPC_TO_CM));
  z /= static_cast<float>(Double(PhysConst::MPC_TO_CM));

  // Add a small epsilon offset (in Mpc) if x or y are exactly zero
  const float epsilon = 1e-6;
  if (x < epsilon) x = epsilon;
  if (y < epsilon) y = epsilon;
  if (z < epsilon) z = epsilon;

  float bx = m_noise_x.GetNoise(x, y, z);
  float by = m_noise_y.GetNoise(x, y, z);
  float bz = m_noise_z.GetNoise(x, y, z);

  Vec3D b_vec(bx, by, bz);
  Vec3D b_field_dir = (b_vec / b_vec.Norm());

  return b_field_dir;
}

} // namespace IGCascade