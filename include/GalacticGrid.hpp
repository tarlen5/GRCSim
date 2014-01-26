/*! \file    GalacticGrid.hpp          \n
             Vec3D class header file
  
  \author    Yusef Shafi	       \n
             UCLA		       \n
	     yshafi@ucla.edu	       \n
  
  \author    Stephen Fegan             \n
             UCLA                      \n
	     sfegan@astro.ucla.edu     \n

  \author    Tim Arlen                 \n
             UCLA                      \n
	     timothyarlen@gmail.com    \n
  
  \date      06/01/2007
  
  \version   1.1
  
  \revision: 1.0 - added code to convert VEC3D_T to double via
             the convert class.
	     1.1 - Add function to do propagation through expanding
	     space.
  \note

*/

#ifndef IGCASCADE_GALACTICGRID_H
#define IGCASCADE_GALACTICGRID_H

#include "Vec3D.hpp"
#include "Vec4D.hpp"
//#include "RandomNumbers.hpp"
#include "RelParticle.hpp"
#include "convert.hpp"
#include "PhysicsConstants.hpp"

#include <iostream>
#include <fstream>
#include<qd/dd_real.h>
//#include<qd/dd.h>
#include <cmath>
#include <map>
#include <cstdlib>
#include <string>

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
//#include <fildes.h>

#include "TRandom3.h"

namespace IGCascade
{

  class ICoord
  {
  public:
    ICoord(int _ix, int _iy, int _iz): ix(_ix), iy(_iy), iz(_iz) { }
    bool operator == (const ICoord& o) { 
      return (ix==o.ix)&&(iy==o.iy)&&(iz==o.iz); 
    }
    bool operator != (const ICoord& o) { 
      return (ix!=o.ix)||(iy!=o.iy)||(iz!=o.iz); 
    }
    int ix;
    int iy;
    int iz;
    
    void DumpShort(std::ostream& stream ) const { 
      stream << ix << ' ' << iy << ' ' << iz; 
    }

  };

  inline bool operator < (const ICoord& c1, const ICoord& c2)
  {
    if(c1.ix < c2.ix)return true;
    else if(c2.ix < c1.ix)return false;
    else if(c1.iy < c2.iy)return true;
    else if(c2.iy < c1.iy)return false;
    else return c1.iz < c2.iz;
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Typedefs///////////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //  typedef Vec3D B;
  typedef std::map<ICoord,Vec3D> Field;

  class MagneticGrid
  {
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Constructors///////////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
    //Default constructor
    inline MagneticGrid();
			
    //Overloaded Constructors
    //MagneticGrid(RandomNumbers * _rng);
    MagneticGrid(TRandom3* _rng, VEC3D_T B_mag, std::string s_cell_size,
		 std::string sfilename);
    //MagneticGrid(RandomNumbers * _rng, VEC3D_T B_mag, std::string s_B_mag, 
    //	 std::string s_cell_size, std::string s_zi);

    //Copy Constructor
    //inline MagneticGrid(const MagneticGrid& _mg);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Member Functions///////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    inline Vec3D CheckMagneticField(const ICoord& c);
    inline Vec3D CheckMagneticField_Lock(const ICoord& c);
    
    ICoord CheckCurrentCell(Vec3D r);

    void PrintFieldMapToFile(std::ofstream& outfile);
    void ReadFieldMap(std::ifstream& infile);
    int LockAttempt(const char* filename, int* fd, std::string rwtype);
    
    void PropagateBFieldRedshift(RelParticle& Photon, RelParticle*&Lepton,
				 Vec3D& n_eo, VEC3D_T& PL, 
				 VEC3D_T& delta_z, const bool LOCK=false);
    
    void Propagation(RelParticle& Photon,RelParticle*& Lepton,VEC3D_T& PL,
		     Vec3D& r_new,VEC3D_T& time_delay,Vec3D& n_eo,
		     VEC3D_T& delta_z,VEC3D_T& delta_zs,Vec3D& e_b);

    Vec3D UpdatePosition(RelParticle*& Lepton,VEC3D_T& PL,Vec3D& n_eo,
			 Vec3D& e_b);

    VEC3D_T Delta_z(VEC3D_T& PL, RelParticle*& Lepton);
    VEC3D_T Delta_zs(VEC3D_T& PL, RelParticle*& Lepton, 
		     Vec3D&r_new, VEC3D_T& delta_time);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Public Data Members////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are now public so we can write them to the output file...
    double         m_cellsize;
    VEC3D_T        m_bmag;  //B_magnitude;
    VEC3D_T        m_DE;     // relative computation precision of roots    
    std::string    m_sfilename;


  private:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Private Data Members////////////////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRandom3* m_rng;
    Field          m_MagneticField;
    
  };


}


namespace std
{

  inline ostream& operator << (ostream& stream, const IGCascade::ICoord& ic);
  inline istream& operator >> (istream& stream, IGCascade::ICoord& ic);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Stream insertion
  inline ostream& operator << (ostream& stream, const IGCascade::ICoord& ic)
  {
    ic.DumpShort(stream);
    return stream;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Stream extraction
  inline istream& operator >> (istream& stream, IGCascade::ICoord& ic)
  {
    char c;
    stream >> c >> ic.ix >> ic.iy >> ic.iz >> c;
    return stream;
  }

}

#endif // IGCASCADE_GALACTICGRID_H
