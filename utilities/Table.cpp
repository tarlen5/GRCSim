/*!
-------------------------------------------------------------------------------
    \file   Table.cpp

    Table class implementation file; classes Table1D and Table2D inherit
    from class Table.

    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
	       arlen@astro.ucla.edu                  \n

    \date      October 10, 2010

    \version:  1.0

    \revision:

    \note:

-------------------------------------------------------------------------------
*/


#include "Table.hpp"


Table1D::Table1D(const std::string& tablefile)
/*!
  Constructor for Table1D.

  \param  tablefile - file for 1D table.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WARNING: required format of tablefile: 1 column of y values,
  followed by 1 column of x values.
  BOTH COLUMNS MUST BE SORTED for future fns to work.
  It is permissible to start a commented line with the hash ('#')
  symbol.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*/
{

  std::ifstream filestream(tablefile.c_str());

  // Check if file exists!!!
  if(!filestream) {
    std::cerr<<"ERROR: "<<__PRETTY_FUNCTION__<<"\n file: "<<tablefile<<
      " cannot be found\n"<<std::endl;
    exit(EXIT_FAILURE);
  }

  const unsigned MAXCHAR = 800;
  char line[MAXCHAR];
  while( !filestream.getline(line,MAXCHAR).eof() ) {
    if (line[0] == '#')  continue;

    std::istringstream ss_line(line);
    double temp_y, temp_x;
    if (filestream) {
      ss_line >> temp_y >> temp_x;
      m_xvector.push_back(temp_x);
      m_yvector.push_back(temp_y);
    }

  }

    // Testing:
    //for(int index=0; index < m_xvector.size(); index++) {
    //std::cout<<"x: "<<m_xvector[index]<<" y: "<<m_yvector[index]<<std::endl;
    //}

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double Table1D::LinInterpolate(const double& xval)
{

  // Perform checks on xval at beginning:
  if(xval < m_xvector[0]) {
    std::cerr<<"\nERROR: "<<__PRETTY_FUNCTION__<<"\nxval cannot be "<<
      "less than (or equal to) the first table value."<<std::endl<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (xval > m_xvector[m_xvector.size()-1]) {
    std::cerr<<"\nERROR: "<<__PRETTY_FUNCTION__<<"\nxval cannot be "<<
      "greater (or equal) to last table value\n"<<std::endl<<std::endl;
    exit(EXIT_FAILURE);
  }


  // Find x0,x1; y0,y1:
  double x0 = 0.0;
  double x1 = 0.0;
  double y0 = 0.0;
  double y1 = 0.0;
  unsigned ix;
  for( ix = 0; ix < m_xvector.size(); ix++) {
    if (m_xvector[ix] > xval) break;
    x0 = m_xvector[ix];
    x1 = m_xvector[ix+1];
    y0 = m_yvector[ix];
    y1 = m_yvector[ix+1];
  }
    
  //std::cout<<"x0: "<<x0<<" x1: "<<x1<<std::endl;
    
  double yval = y0 + (xval - x0)*(y1 - y0)/(x1 - x0);
  //std::cout<<"yval: "<<yval<<std::endl;
    
  return yval;
}



/*********************************************************\
 * class Table2D                                         *
 *                                                       *
 *                                                       *
\*********************************************************/

Table2D::Table2D(const std::string& tablefile)
/*!
  Constructor for Table2D.

  \param  tablefile - file for 2D table.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WARNING: required format of tablefile: for f(x1,x2)
  need n1 + 1 rows of variable x1, and n2 + 1 columns of x2
  where n1 is number of x1 points on x1 axis & n2 is num of
  x2 points on x2 axis.

  NOTE: first line MUST be -1 x2_1 x2_2 ...
  all subsequent lines must be x1_i f(x1_i,x2_j) ...

  BOTH AXES MUST BE SORTED for table class fns to work.
  It is permissible to start a commented line with the hash ('#')
  symbol.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*/
{

  std::ifstream filestream(tablefile.c_str());
  // Check if file exists.
  if(!filestream) {
    std::cerr<<"ERROR: "<<__PRETTY_FUNCTION__<<"\n file: "<<tablefile<<
      " cannot be found\n"<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  int irow = 0; int icol = 0;
  while(getline(filestream,line)) {
    if (line[0] == '#')  continue;
    m_table.push_back( std::vector<double>() );
    std::istringstream ss_line(line);
    icol = 0;
    while(ss_line) {
      double temp_table_val;
      ss_line >> temp_table_val;
      if(ss_line) {
	//std::cout<<temp_table_val<<" ";
	m_table[irow].push_back(temp_table_val);
	icol++;
      }
      //std::cout<<"\n";
    }
    irow ++;

    // Checking format of file:
    if(irow==1) {
      m_col = icol;
    } else {
      if(m_col != icol) {
	std::cerr<<"\nERROR: in "<<__PRETTY_FUNCTION__<<" incorrect format"<<
	  " of file: "<<tablefile<<". Number of columns must be the same"<<
	  " for all rows.\n"<<std::endl;
      }
    }
    //////////////////////////////

  }
  m_row = irow;

  //DisplayTable();
  //std::cout<<"rows: "<<m_row<<" cols: "<<m_col<<std::endl;

}


Table2D::Table2D(std::vector<double>& rows_vec,
		 std::vector<double>& cols_vec)
{

  //std::cout<<"\n  cols: "<<cols_vec.size()<<std::endl;
  //std::cout<<"  rows: "<<rows_vec.size()<<std::endl;

  ////////////////////////////////
  /// DO ALL COLS ...
  // Do first col separately:
  m_table.push_back( std::vector<double>() );
  m_table[0].push_back(-1);   // [0][0] element is -1
  for(unsigned icol=0; icol<cols_vec.size(); icol++)
    m_table[0].push_back(cols_vec[icol]);
  m_col = m_table[0].size();
  //std::cout<<"m_col: "<<m_col<<std::endl;
  // First row complete...

  // Now implement all rows:
  for(unsigned irow=0; irow<rows_vec.size(); irow++) {
    m_table.push_back( std::vector<double>() );
    m_table[irow+1].push_back(rows_vec[irow]);
    for(unsigned icol=1; icol<m_col; icol++) {
      m_table[irow+1].push_back(0.0);
    }
  }
  m_row = m_table.size();
  //std::cout<<"m_row: "<<m_row<<std::endl;
  //DisplayTable();

}


void Table2D::DisplayTable(void)
{
  // Output table:
  for(int irow=0; irow<m_row; irow++) {
    for(int icol=0; icol<m_col; icol++) {
      std::cout<<m_table[irow][icol]<<" ";
    }
    std::cout<<std::endl;
  }

}
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double Table2D::LinInterpolate(const double& x1val, const double& x2val)
/*!
  Performs bilinear interpolation on the grid square. If you need a refresher
  on what this is, see Numerical Recipes for C, section 3.6.

  \param - x1val Desired x1 coordinate f(x1,x2) (row in table)

  \param - x2val Desired x2 coordinate f(x1,x2) (col in table)

  WARNING: does not check that table is ordered as low -> high across cols
  and low -> high across rows.

*/
{

  ////////////////////////////////////////////////////////////////////////
  // ERROR checking. Make sure we're in range...
  // Again, WARNING: m_table must be ordered as low -> high across row
  //   and low -> high down column.
  ////////////////////////////////////////////////////////////////////////
  double x1_min = m_table[1][0];
  double x1_max = m_table[m_row-1][0];
  if ( (x1val < x1_min) || (x1val > x1_max) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "x1val not in range.\n" <<
      "x1val: "<<x1val<<std::endl <<
      "x1_min: "<<x1_min<<std::endl <<
      "x1_max: "<<x1_max<<std::endl;
    exit(EXIT_FAILURE);
  }

  double x2_min = m_table[0][1];
  double x2_max = m_table[0][m_col-1];
  if( (x2val < x2_min) || (x2val > x2_max) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "x2val not in range.\n" <<
      "x2val: "<<x2val<<std::endl <<
      "x2_min: "<<x2_min<<std::endl <<
      "x2_max: "<<x2_max<<std::endl;
    exit(EXIT_FAILURE);
  }
  ////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  // Define points 1-4
  // Assumes x1val and x2val fall between the points on the grid defined by
  //   y1=m_table[j][k], y2=m_table[j+1][k], y3=m_table[j+1][k+1],
  //   y4m_table[j][k+1]
  //////////////////////////////////////////////////////////////////////////
  int col_lo = 0; int row_lo=0;
  double x1lo=0.0; double x1hi=0.0;double x2lo=0.0; double x2hi=0.0;
  for(int irow = 1; irow < m_row; irow++) {
    if(m_table[irow][0] < x1val) continue;
    else {
      x1hi = m_table[irow][0];
      x1lo = m_table[irow-1][0];
      row_lo = irow-1;
      break;
    }
  }
  //std::cout<<"x1lo: "<<x1lo<<" x1hi: "<<x1hi<<std::endl;
  for(int icol = 1; icol < m_col; icol++) {
    if(m_table[0][icol] < x2val) continue;
    else {
      x2hi = m_table[0][icol];
      x2lo = m_table[0][icol-1];
      col_lo = icol-1;
      break;
    }
  }
  //std::cout<<"x2lo: "<<x2lo<<" x2hi: "<<x2hi<<std::endl;


  // Do interpolation:
  double interp_t = (x1val-x1lo)/(x1hi - x1lo);
  double interp_u = (x2val-x2lo)/(x2hi - x2lo);

  double y1 = m_table[row_lo][col_lo];
  double y2 = m_table[row_lo+1][col_lo];
  double y3 = m_table[row_lo+1][col_lo+1];
  double y4 = m_table[row_lo][col_lo+1];
  //std::cout<<"y1: "<<y1<<" y2: "<<y2<<" y3: "<<y3<<" y4: "<<y4<<std::endl;
  double val = (1-interp_t)*(1-interp_u)*y1 + interp_t*(1-interp_u)*y2 +
    interp_t*interp_u*y3 + (1-interp_t)*interp_u*y4;

  return val;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double Table2D::GetVal(const int row, const int col)
{

  if ( (row < 0) || (row >= m_row) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "row not in range!\n" <<
      "  row: "<<row<<std::endl <<
      "  row_max: "<<m_row<<std::endl;
    exit(EXIT_FAILURE);
  }

  if ( (col < 0) || (col >= m_col) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "col not in range!\n" <<
      "  col: "<<col<<std::endl <<
      "  col_max: "<<m_col<<std::endl;
    exit(EXIT_FAILURE);
  }

  return m_table[row][col];

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double Table2D::GetRowVal(const int row)
{
  if ( (row < 0) || (row >= m_row) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "row not in range!\n" <<
      "  row: "<<row<<std::endl <<
      "  row_max: "<<m_row<<std::endl;
    exit(EXIT_FAILURE);
  }

  return m_table[row][0];
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double Table2D::GetColVal(const int col)
{

  if ( (col < 0) || (col >= m_col) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "col not in range!\n" <<
      "  col: "<<col<<std::endl <<
      "  col_max: "<<m_col<<std::endl;
    exit(EXIT_FAILURE);
  }


  return m_table[0][col];
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Table2D::SetVal(const int row, const int col, const double val)
{
  if ( (row < 1) || (row >= m_row) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "row not in range!\n" <<
      "  row: "<<row<<std::endl <<
      "  row_max: "<<m_row<<std::endl;
    exit(EXIT_FAILURE);
  }
  if ( (col < 1) || (col >= m_col) ) {
    std::cerr << "ERROR: in function: "<<__PRETTY_FUNCTION__<<"\n";
    std::cerr << "col not in range!\n" <<
      "  col: "<<col<<std::endl <<
      "  col_max: "<<m_col<<std::endl;
    exit(EXIT_FAILURE);
  }


  m_table[row][col] = val;

}

