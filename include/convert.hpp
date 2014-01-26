/*! \file   convert.hpp
            Converts type VEC3D_T (a quad double or double double
	    precision number) to type double.
  
  \author   Tim Arlen              \n
            UCLA                   \n
	    timothyarlen@gmail.com \n

  \date     5/02/2008
  \version  1.0
            1.1 - Name changed from convertToDouble() to Double()
  \note

*/

#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <qd/dd_real.h>

typedef dd_real VEC3D_T;

// NOTE: Should I add this code to the qd/dd_real class?

namespace IGCascade
{
  
  class BadConversion : public std::runtime_error {
  public:
    BadConversion(const std::string& s)
      : std::runtime_error(s)
    { }
  };


  inline double Double(VEC3D_T x) {
    std::ostringstream o;
    if (!(o<<x))
      throw BadConversion("Double(double)");

    std::string s = o.str();
    std::istringstream i(s);
    double y;

    if (!(i>>y))
      throw BadConversion("Double(\"" + s + "\")");
    
    return y;

  }
  
}

#endif // ifndef CONVERT_HPP
