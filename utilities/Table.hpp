/*! 
-------------------------------------------------------------------------------
  \file    Table.hpp
           Table class header file
  
  \author   Timothy C. Arlen       \n
            UCLA                   \n
	    arlen@astro.ucla.edu   \n

  \date     10/08/10
  \version  0.0
  \note     Works only with double, not dd_real.
------------------------------------------------------------------------------  
*/


#ifndef IGCASCADE_TABLE_H
#define IGCASCADE_TABLE_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <iomanip>

//#include "PhysicsConstants.hpp"
//#include "convert.hpp"

  
/////////////////////////////////////////////////
/// Table class: Abstract Base Class
/////////////////////////////////////////////////
class Table
{
public:
  Table() { /* nothing */ }
  virtual ~Table() { /*nothing*/ }
  
protected:
  
  //virtual void LinInterpolate();
  virtual void DisplayTable(void) = 0;
  
};


class Table1D : public Table
{
public:
  
  /////// Constructor /////////
  Table1D(const std::string& tablefile);
  
  /////// Functions /////////
  double LinInterpolate(const double& x1val);
  
private:
  
  /////// Functions /////////
  void DisplayTable(void) { /* nothing yet but define later*/ }
  
  /////// Variables ////////
  std::vector<double> m_xvector;
  std::vector<double> m_yvector;    
  
  
};

class Table2D: public Table
/*
  Possible improvements:
  - cubic spline interpolation
  - 
*/
{
public:
  
  ////// Constructor //////
  Table2D(const std::string& tablefile);
  Table2D(std::vector<double>& rows_vec, 
	  std::vector<double>& cols_vec);
  
  ////// Functions //////
  double LinInterpolate(const double& x1val, const double& x2val);
  void   DisplayTable(void);

  //////////// Accessors /////////////
  int    GetNRows() { return m_row; }
  int    GetNCols() { return m_col; }
  double GetRowVal(const int row);
  double GetColVal(const int col);
  double GetVal(const int row, const int col);

  //////////// Simple Settors /////////
  void SetVal(const int row, const int col, const double val);

  std::vector< std::vector<double> > m_table;
  
private:
  
  ////// Functions //////  
  int m_col;                   // Num cols in table (including the axes)
  int m_row;                   // Num rows in table
  
};


#endif // IGCASCADE_TABLE_H
