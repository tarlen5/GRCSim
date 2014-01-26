//We will move this file's functionality to KleinNishina class eventually

//-*-mode:c++; mode:font-lock;-*-

/*! \file Vec3D.hpp
  Vec3D class header file
  
  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n
  \author	Yusef Shafi
			UCLA
		yshafi@ucla.edu

  \date     1/24/2007
  \version  0.2
  \note
		
*/


#ifndef PHYSICS_CONSTANTMAGNETICFIELD_H
#define PHYSICS_CONSTANTMAGNETICFIELD_H

#include <map>
#include<qd/dd_real.h>
//#include<qd/dd.h>

#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "RandomNumbers.hpp"

namespace Physics
{
  /*
    class ICoord
    {
    public:
    ICoord(int _ix, int _iy, int _iz): ix(_ix), iy(_iy), iz(_iz) { }
    int ix;
    int iy;
    int iz;
    };

    inline operator < (const ICoord& c1, const ICoord& c2)
    {
    if(c1.ix < c2.ix)return true;
    else if(c2.ix < c1.ix)return false;
    else if(c1.iy < c2.iy)return true;
    else if(c2.iy < c1.iy)return false;
    else return c1.iz < c2.iz;
    }

    typedef Vec3D B;
    typedef map<ICoord,B> Field;

    Vec3D CheckMagneticField(Vec3D &r, Field &bfield, RandomNumbers &rng)
    {

 
    double x = r.x;
    double y = r.y;
    double z = r.z;
    ICoord c(x, y, z);
    
    if(bfield.find(c) == bfield.end())
    {
    // cell has not been used before so create the field
    //B bnew = generateRandomField();	  
    B bnew;
    bnew = bnew.UniformSphereDirection(&rng);
    bfield[c] = bnew;
    }
      
    B b = bfield[c];
	
    return b;
    }*/

}

#endif
