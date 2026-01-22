// Test for MFTurbulentContinuous: Propagate a mock lepton through a turbulent
// magnetic field and log the B-field direction at each step to a CSV file.
// Usage: ./test_MFTurbulentContinuous <B_field_G> <coh_len_Mpc> <energy_eV>
// <z_init>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "HighPrecProp.hpp"
#include "MFTurbulentContinuous.h"
#include "PhysicsConstants.hpp"
#include "RelParticle.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"

using namespace IGCascade;

float rv_propagation_length(RandomNumbers *rng, float distance_mpc) {
  // Uniform distribution between 0 and coh_len_Mpc / 10
  float rand_uniform = static_cast<float>(rng->Uniform());
  return rand_uniform * distance_mpc;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " <B_field_G> <coh_len_Mpc> <energy_GeV> <z_init>\n";
    return 1;
  }
  double B_field_G = std::stod(argv[1]);
  double coh_len_Mpc = std::stod(argv[2]);
  double energy_GeV = std::stod(argv[3]);
  double z_init = std::stod(argv[4]);

  // Set up random number generator for propagation distance
  RandomNumbers *rng = new RandomNumbers(0.0, 1.0);
  auto random_dist = [coh_len_Mpc, rng]() {
    float rand_uniform = static_cast<float>(rng->Uniform());
    return rand_uniform * static_cast<float>(coh_len_Mpc / 10.0);
    // return rv_propagation_length(rng, static_cast<float>(coh_len_Mpc / 10.0));
  };

  // Set up RelParticle electron as a pointer
  RelParticle *lepton = new RelParticle();
  lepton->m_r4 = Vec4D(0.0, 0.0, 0.0, 0.0);
  VEC3D_T E_e = energy_GeV * 1.0E9; // Convert GeV to eV
  Vec3D n1(0.0, 0.0, 1.0); // initial direction, will be randomized each step
  lepton->m_p4.r0 = E_e;
  lepton->m_m0 = PhysConst::SI_MELEC;
  lepton->m_p4.r = n1 * E_e;
  lepton->m_z = z_init;
  lepton->m_z_s = z_init;
  lepton->m_z_s_int = z_init;
  lepton->m_z_int = z_init;
  lepton->m_rnext = 0.0;
  lepton->m_q = -1;
  lepton->m_weight = 1.0;

  // Set up dummy photon (only required values for PropagateConstantMF)
  RelParticle photon;
  photon.m_r4 = lepton->m_r4;
  photon.m_p4.r0 = E_e; // arbitrary, just needs to be nonzero
  photon.m_p4.r = n1 * E_e;
  photon.m_z = z_init;
  photon.m_z_s = z_init;

  // Set up magnetic field
  VEC3D_T b_mag = B_field_G;
  VEC3D_T coh_len = coh_len_Mpc;
  MFTurbulentContinuous mf(rng, b_mag, coh_len);

  // Open CSV file
  // Build output filename with parameters
  std::string out_filename =
      "tests/MFTurbulence_PropLength_vs_BDirection_" + std::string(argv[1]) +
      "G_" + std::string(argv[2]) + "Mpc_" + std::string(argv[3]) + "GeV_" +
      std::string(argv[4]) + "z.csv";

  std::ofstream fout(out_filename);
  fout << "Prop_length_Mpc,pos_x,pos_y,pos_z,e_b_x,e_b_y,e_b_z\n";

  VEC3D_T propagation_length_taken = 0.0;
  unsigned int nstep = 0;
  unsigned int max_steps = 1000;
  // Main propagation loop
  while (lepton->m_z_s > 0.0) {
    // Get B-field direction at current position
    Vec3D b_dir = mf.GetMagneticFieldDirectionAtPosition(lepton->m_r4.r);

    // Step: choose random propagation length (in Mpc, convert to cm)
    double prop_len_Mpc = random_dist();
    double prop_len_cm = prop_len_Mpc * Double(PhysConst::MPC_TO_CM);
    VEC3D_T prop_len_cm_hp = prop_len_cm; // Convert to high-precision type

    // Generate a random direction for this step
    Vec3D n_step = Vec3D::UniformSphereDirection(rng);

    // Calculate new position
    Vec3D r_new = lepton->m_r4.r + n_step * prop_len_cm;

    // Dummy time delay and redshift changes
    VEC3D_T delta_time = 0.0;
    VEC3D_T delta_z = 0.0;
    VEC3D_T delta_zs = 0.0;

    // Use PropagateConstantMF to update lepton and photon
    PropagateConstantMF(
        photon, lepton, b_mag, prop_len_cm_hp, r_new, delta_time, n_step,
        delta_z, delta_zs, b_dir
    );

    propagation_length_taken += prop_len_Mpc;
    fout << propagation_length_taken << "," << lepton->m_r4.r.x << ","
         << lepton->m_r4.r.y << "," << lepton->m_r4.r.z << "," << b_dir.x << ","
         << b_dir.y << "," << b_dir.z << "\n";

    // For this first test, just take 100 steps:
    nstep++;
    if (nstep >= max_steps) {
      break;
    }
  }

  fout.close();
  std::cout << "Test complete. B-field directions logged to " << out_filename
            << "\n";

  delete lepton;
  delete rng;

  return 0;
}
