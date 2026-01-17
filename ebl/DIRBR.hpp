/*! \file    DIRBR.hpp

    Diffuse InfraRed Background Radiation

    \author Vladimir V. Vassiliev \n
            Department of Physics and Astronomy \n
            UCLA \n
            E-Mail: vvv@astro.ucla.edu

    \date   November 13, 2004

    \version 1.0

    \note

 */

/*! \class  DIRBR
    \brief  Diffuse InfraRed Background Radiation

*/

#ifndef DIRBR_H
#define DIRBR_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
// #include <cmath>
#include "DIRBRBase.hpp"

class DIRBR : public DIRBRBase {
public:
  DIRBR();               //!< constructor
  DIRBR(double, double); //!< overloaded constructor
  virtual ~DIRBR();      //!< destructor

  //! Get... accessors
  //@{
  virtual double GetDIRBR(double lambda, double z = 0.0); //!< DIRBR
  virtual double GetCMBR(double &lambda);                 //!< CMBR
  double GetEnergy(double, double); //!< Energy within wavelength interval
  double GetDIRBRPowerLaw(double);  //!< Power Law DIRBR
  //@}

  //! Set... accessors
  //@{
  DIRBR_ERR SetDIRBR(int, double *, double *); //!< set DIRBR model
  void SetDIRBRPowerLaw(double);               //!< set power law DIRBR
  //@}

  //! methods
  //@{
  DIRBR_ERR TestDIRBRlimits(void);             //!< test existing direct limits
  double OpticalDepthPowerLaw(double, double); //!< find opt depth
  virtual double OpticalDepth(double &E,
                              double &z); //!< find opt depth for E, z
  virtual double Distance(double &z);     //!< find distance to redshift z
  //@}

protected:
  //!
  //@{

  //@}

private:
  //!
  //@{
  int K;            //!< number of specified wavelength bands
  double *I;        //!< pointer to SED array
  double *lambda;   //!< pointer to lambda array
  double *LnI;      //!< pointer to Ln of SED array
  double *Lnlambda; //!< pointer to Ln of lambda array
  double *h;        //!< pointer to array of spline intervals
  double *m;        //!< pointer to spline array
  double dlnlambda; //!< ln(wavelength) integration step

  double sigma240; //!< sigma for FIRAS 240 \mum SED
  double Imax240;  //!< MAX SED at FIRAS 240 \mum
  double Imin240;  //!< MIN SED at FIRAS 240 \mum

  double Imax1216A; //!< MAX SED at 1216 A
  double Imin1216A; //!< MIN SED at 1216 A
  double I1216A;    //!< SED at 1216 A
  double IndexUV;   //!< SED spect ind in extreme UV region 900A<...<1216A

  double nuFnu; //!< SED scaling constant [nW m^-2 sr^-1]
  double kappa; //!< wavelength conversion constant

  double q_min; //!< range of definition for optical depth integration
  int Nq_min;   //!< number of points in intrinsic kernel arrays
  double *q;    //!< auxiliary array for calculating optical depth
  double *f;    //!< auxiliary array for calculating optical depth
  double *F3;   //!< auxiliary array for calculating optical depth

  double dz; //!< redshift integration step;

  double gamma;    //!< spectral index for power law SED
  double PLfactor; //!< SED factor for power law DIRBR

  DIRBR_ERR DIRBR_err; //!< DIRBR error code
                       //@}

  //! methods
  //@{
  DIRBR(const DIRBR &);                  //!< blocking copy constructor
  const DIRBR &operator=(const DIRBR &); //!< blocking copy constructor
  void AllocateMemoryT();        //!< Allocates and calculates internal arrays
  void DeleteMemoryT();          //!< Removes internal tau arrays
  void DeleteMemoryS();          //!< Removes internal spline arrays
  double SEDFactor(void);        //!< Power law specific e and theta integral
  double RedShiftFactor(double); //!< Redshift factor
  //@}

  //! Get.. accessors
  //@{
  double GetFIRAS(double); //!< COBE FIRAS result
  //@}

  //! Set.. accessors
  //@{
  DIRBR_ERR SetFIRAS240(double); //!< Sets SED at 240 \mum
  DIRBR_ERR Set1216A(double);    //!< Sets SED at 1216 A
  //@}
};

#endif // DIRBR_H
