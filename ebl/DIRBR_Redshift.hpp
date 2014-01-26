/*! 

  \file DIRBR_Redshift.hpp
        DIRBR_Redshift class header file
  
  \author   Tim Arlen               \n
            UCLA                    \n
            arlen@physics.ucla.edu  \n
	    
  \date     July 24, 2010
  \version  0.0
  \note

*/

#ifndef DIRBR_REDSHIFT_H
#define DIRBR_REDSHIFT_H

#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<vector>

//#include "PhysicsConstants.hpp"
//#include "convert.hpp"
#include "DIRBRBase.hpp"
#include "Table.hpp"

class DIRBR_Redshift: public DIRBRBase
{
public:
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Constructors///////////////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DIRBR_Redshift(std::string filename);
  

  virtual ~DIRBR_Redshift() {
    /* nothing to see... */
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Public Member Functions////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  virtual double GetDIRBR( double lambda, double z);
  virtual double GetCMBR( double& lambda);
  virtual double OpticalDepth( double& E, double& z );
  virtual double Distance( double& z ); 

  virtual DIRBR_ERR SetDIRBR(int,double*,double*);     //!< set DIRBR model 
  virtual DIRBR_ERR TestDIRBRlimits(void); //!< test existing direct limits

private:

  void InitializeIntegral(void);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Private Data Members //////////////////////////////
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Table2D* SEDTable;
  
  double m_lambda_min;
  double m_lambda_max;
  double m_zmin;
  double m_zmax;
  
  /////////////////////////////////////
  // Integration parameters
  /////////////////////////////////////
  int m_nsteps;
  double m_qmin;
  double m_grid_factor;
  double m_dlnq;
  double* m_Farray;
  double* m_qarray;
  
  double m_dz;

};


#endif  // DIRBR_REDSHIFT_H

