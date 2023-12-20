/*! \file KleinNishina.hpp
  KleinNishina class header file
  
  
  \author   Yusef Shafi            \n
            UCLA                   \n
	    yshafi@ucla.edu        \n

  \author   Tim Arlen              \n
            UCLA                   \n
	    timothyarlen@gmail.com \n

  \author   Vladimir Vassiliev     \n
            UCLA                   \n
	    vvv@astro.ucla.edu     \n

  \date     08/17/2006
  \version  0.0
  \note
  Think of this class as a "space" where +
  various scattering processes "operate"
  on a photon and electron four vector.
*/

#ifndef IGCASCADE_KLEINNISHINA_H
#define IGCASCADE_KLEINNISHINA_H

#include<vector>
#include<string>
#include<cmath>

#include<qd/dd_real.h>

#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "RandomNumbers.hpp"
//#include "RelParticle.hpp"
#include "DIRBRBase.hpp"
#include "convert.hpp"
#include "PhysicsConstants.hpp"

// For Star Structure:
typedef struct {
	VEC3D_T Radius;                       // [cm]
	VEC3D_T Temp;                         // [eV]
	VEC3D_T Luminosity;                   // [eV/s]
} StructStar;


namespace IGCascade
{

  class KleinNishina
  {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Constructors///////////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
    //Overloaded Constructor
    /// \param: RandomNumbers reference
    KleinNishina(RandomNumbers* _rng);
			
    //Copy Constructor... are we sure we are checking for self assignment?
    //! note: default constructor is provided as a
    //	formality but should probably NOT be used
    /// \param: KleinNishina reference
    inline KleinNishina(const KleinNishina& _kn);
			
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Static Functions///////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
    //Maybe implement scattering processes as static methods

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Member Functions///////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    //TO BE MADE PRIVATE (only here for testing)
    VEC3D_T RelativisticKinematics(Vec4D& P_e, Vec4D& P_p);
    VEC3D_T PropagationLengthBB(VEC3D_T Ebb, Vec4D& P_e, Vec4D& P_p);
    VEC3D_T PropagationLengthCMBandIR(DIRBRBase* ebl_model, VEC3D_T Ebb, 
				      Vec4D& P_e, Vec4D& P_p);
    bool PropagationLengthStar(StructStar& star, Vec4D& P_e, 
			       Vec4D& R_e, Vec4D& P_p, 
			       Vec4D& R_p, VEC3D_T& pl);
    //Dump Function to get info out
    //void Dump(std::ostream& stream = std::cout) const;
			
    			
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Public Data Members////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		
		
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Private Member Functions///////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
    //Default constructor
    /// \param: N/A
    inline KleinNishina();
			
    //To be called by ALL scattering methods
						
    VEC3D_T Scattering(VEC3D_T x);
    VEC3D_T IncompleteTotalCrossSection(VEC3D_T z, VEC3D_T x);
    VEC3D_T IsotropicRadiation(VEC3D_T  z, VEC3D_T gamma);
    VEC3D_T IsotropicSigma(VEC3D_T  z, VEC3D_T gamma);
    void    IsotropicFactors(VEC3D_T & f1, VEC3D_T & f2, VEC3D_T x);


    //VEC3D_T natLog(VEC3D_T x);
    VEC3D_T PolyLog1(VEC3D_T x);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Private Data Members///////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  
    RandomNumbers* m_rng;

    // Relevant mathematical constants
    VEC3D_T m_PI2d6;
	  
    // computation accuracy
    VEC3D_T m_DE;
	  
    // Isotropic Radiation section 
    // Saves time for multiple scattering of 
    // photons with the same energy 
    VEC3D_T m_IRd_z;
    VEC3D_T m_IRd_b;
    VEC3D_T m_IRd_f1c;
    VEC3D_T m_IRd_f2c;
    VEC3D_T m_IRd_f1b;
    VEC3D_T m_IRd_f2b;
    VEC3D_T m_IRd_f1a;
    VEC3D_T m_IRd_f2a;
    VEC3D_T m_IRd_U;

    // BB Propagation length section
    // Parameters for numerical integration
    VEC3D_T m_egy_bb;        // kT of BB spectrum [eV]
    VEC3D_T m_num_dens_bb;   // photon number density [cm^-3]
    VEC3D_T m_egy_dens_bb;   // photon energy density [eV cm^-3]

    double m_dx;  // integration step from 0 to xU
    //double m_xU;  // energy upper bound; physical units are xU*kT
    int m_num_int;    // number of integration steps (int)((double) xU/dx)

    // For Integrals in PropagationLengthCMBandIR():
    VEC3D_T egy_min;   // [eV]
    VEC3D_T egy_max;     // [eV]
    VEC3D_T m_num_bins_dec;
    VEC3D_T m_egy_factor;
    std::vector<double> m_lambda_vec;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Default class constructor
  inline KleinNishina::KleinNishina()
  {
    m_rng = new RandomNumbers(0.0, 1.0);  //("random_numbers.seed");
  }
	
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Overloaded class constructor
  /// \param _rng: rng
  /*inline KleinNishina::KleinNishina(RandomNumbers * _rng): rng(_rng)
    {
    // nothing to see here
    }*/

	
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Copy constructor
  // \param _kn: KleinNishina Object to copy
  inline KleinNishina::KleinNishina(const KleinNishina& _kn): m_rng(_kn.m_rng)
  {
    // nothing to see here
  }

}
#endif // PHYSICS_KLEINNISHINA_H
