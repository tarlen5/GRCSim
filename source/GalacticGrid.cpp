/*!
-------------------------------------------------------------------------------
  \file     GalacticGrid.cpp
            Propagates charged leptons and photons
            through cosmological magnetic fields.

  \author   Timothy C. Arlen                     \n
            Department of Physics and Astronomy  \n
            UCLA                                 \n
            arlen@astro.ucla.edu                 \n

  \author   Yusef Shafi                          \n
            UCLA	                               \n
            yshafi@ucla.edu	                     \n
  \author   Stephen Fegan                        \n
            UCLA                                 \n
            sfegan@astro.ucla.edu                \n

  \date   03/14/2008

  \version 1.1

  \revision
      04/24/2008 - Added Magnetic Field Propagation in Expanding universe
                   in PropagateMagFieldExpansion() function.
      02/06/2009 - Completely changed the way I determined if
                   propagation would take electron to a new cell.

  \note
      1) IMPORTANT: This function only works properly for
        - cell_size >= 0.1 Mpc, and for
        - B >= 10E-10
        - Electron Energy > 0.1 TeV
      due to issues in correctly computing cell-crossing for electron IC
      propagation steps in the cmb photon field.
-------------------------------------------------------------------------------
*/

#include "GalacticGrid.hpp"

using namespace PhysConst;

namespace IGCascade {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Overloaded class constructor (1)
/// \param _MagneticField: m_MagneticField
/// \param _rng: rng
// MagneticGrid(
//     RandomNumbers *_rng, VEC3D_T B_mag, std::string s_cell_size,
//     std::string sfilename
// );

MagneticGrid::MagneticGrid(
    RandomNumbers *_rng, const std::string &B_mag,
    const std::string &s_cell_size, const std::string &redshift,
    const std::string &mf_dir
) {

  // public member:
  m_DE = "1.0E-25";
  m_rng = _rng;

  // cellsize defined in units of Mpc but converted to cm
  std::istringstream(s_cell_size) >> m_cellsize; // [Mpc]
  // MpcToCm = 3.086E+24;                        // [cm/Mpc]
  m_cellsize = m_cellsize * Double(MPC_TO_CM); // [cm]
  m_bmag = B_mag.c_str();                      // [gauss]

  m_sfilename = DefineMFfile(mf_dir, B_mag, s_cell_size, redshift);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Private Methods:
///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::string MagneticGrid::DefineMFfile(
    const std::string &mf_dir, const std::string &s_Bmag,
    const std::string &s_cellsize, const std::string &redshift
) {
  std::string MFfilename = mf_dir + "MagneticGrid_B" + s_Bmag + "_L" +
                           s_cellsize + "_z" + redshift + ".txt";

  // Ensure the directory exists
  std::string dir = mf_dir;
  if (!dir.empty() && dir.back() == '/') dir.pop_back();
  system(("mkdir -p " + dir).c_str());

  // Create the file if it does not exist
  std::ifstream check(MFfilename.c_str());
  if (!check.good()) {
    std::ofstream create(MFfilename.c_str());
    create.close();
    std::cout << "Creating magnetic field file: " << MFfilename << std::endl;
  }
  return MFfilename;
}

void MagneticGrid::PropagateBFieldRedshift(
    RelParticle &Photon, RelParticle *&Lepton, Vec3D &n_eo, VEC3D_T &PL,
    VEC3D_T &delta_z, const bool LOCK
)
/* Designed to work with the KleinNishina class. A propagation
   length has been determined and the IC scattered photon and
   lepton are computed, then the lepton and photon are translated
   and rotated due to the presence of the cosmological magnetic field.

   \param:

     Photon - The photon which will be rotated and translated

     Lepton - The lepton which will be rotated and translated

     n_eo - Lepton direction before IC scat with photon.

     PL - Propagation length of lepton                     [cm]

     delta_z - change in redshift of the lepton over the
               distance of the propagation length.

     LOCK - if 'true', use CheckMagneticField_Lock; default value is 'false'.

     IMPORTANT: assumes delta_z is always a positive quantity.

    \note: returns void, but updates photon/lepton 4-position and
           momentum directly
  */
{

  VEC3D_T PL_original = PL;
  VEC3D_T MIN_Precision = 1.E-12 * PL_original;
  VEC3D_T PL_remain = PL;
  while (PL_remain > MIN_Precision) {

    PL_remain = PL;

    ICoord c_cur = CheckCurrentCell(Lepton->m_r4.r);

    Vec3D e_b(0., 0., 0.);
    if (LOCK)
      e_b = CheckMagneticField_Lock(c_cur);
    else
      e_b = CheckMagneticField(c_cur);

    Vec3D r_new = UpdatePosition(Lepton, PL, n_eo, e_b);
    ICoord c_new = CheckCurrentCell(r_new);

    if (c_new != c_cur) { // Cell crossed
      // std::cerr<<"Cell Crossing"<<std::endl;
      //  Find Length to new cell:
      VEC3D_T PL_left = 0;
      VEC3D_T PL_right = PL;
      while ((PL_right - PL_left) > MIN_Precision) {
        VEC3D_T PL_test = 0.5 * (PL_left + PL_right);
        r_new = UpdatePosition(Lepton, PL_test, n_eo, e_b);
        ICoord c_test = CheckCurrentCell(r_new);
        if (c_test == c_cur)
          PL_left = PL_test;
        else
          PL_right = PL_test;
      }
      // VEC3D_T PL_fraction = PL_right/PL;
      PL = PL_right;
    }

    VEC3D_T delta_z = GetLeptonDelta_z(PL, Lepton);
    VEC3D_T delta_time = 0.0;
    VEC3D_T delta_zs = GetPhotonDelta_zs(PL, Lepton, r_new, delta_time);

    if ((Lepton->m_z_s - delta_zs) >= m_DE) { // z_s = 0 surface not crossed

      PropagateConstantMF(
          Photon, Lepton, m_bmag, PL, r_new, delta_time, n_eo, delta_z,
          delta_zs, e_b
      );
      PL_remain -= PL;

    } else { // z_s = 0 surface crossed

      /////////////////////////////////////////////
      // compute PL to z_s = 0 surface:          //
      /////////////////////////////////////////////
      VEC3D_T PL_right = PL;
      VEC3D_T PL_left = "0.0";
      Vec3D r_e = Lepton->m_r4.r;
      // Compute PL to z_s = 0.0 surface.
      while (fabs(Lepton->m_z_s - delta_zs) > m_DE) // About 40 iterations.
      {
        PL = (PL_right + PL_left) / 2.0;
        r_new = UpdatePosition(Lepton, PL, n_eo, e_b);
        delta_zs = GetPhotonDelta_zs(PL, Lepton, r_new, delta_time);

        if ((Lepton->m_z_s - delta_zs) > 0.0)
          PL_left = PL;
        else
          PL_right = PL;
      }
      /////////////////////////////////////////////

      delta_z = GetLeptonDelta_z(PL, Lepton);

      PropagateConstantMF(
          Photon, Lepton, m_bmag, PL, r_new, delta_time, n_eo, delta_z,
          delta_zs, e_b
      );
      PL_remain = "0.0";
    }

  } // end while
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vec3D MagneticGrid::UpdatePosition(
    RelParticle *&Lepton, VEC3D_T &PL, Vec3D &n_eo, Vec3D &e_b
)
/*!
  Computes the new position of Lepton if propagated a distance PL.
    IMPORTANT: does not change any of Lepton's values or PL.

  \param PL     - prop length of lepton [cm]
  \param Lepton - charged lepton which is to be propagated

  \return r_new - new position of Lepton
*/

{

  VEC3D_T ze = Lepton->m_z;
  // VEC3D_T m_eV = sqrt(Lepton->m_p4*Lepton->m_p4);
  // std::cerr<<"me_eV (Update): "<<m_eV<<std::endl;
  VEC3D_T p_o = Lepton->m_p4.r.Norm() / (1.0 + ze);
  // VEC3D_T rho  = m_eV*m_eV/p_o/p_o;              // [unitless]

  VEC3D_T Kappa_o = CGS_C * 1.E-8;
  Kappa_o *= Lepton->m_q * m_bmag / p_o;

  VEC3D_T KappaL = Kappa_o * PL;

  Vec3D A = n_eo * PL;
  Vec3D DeltaR = A;
  // Play with m_bmag values to see what best precision is:
  if (fabs(KappaL) > 0.0) {
    Vec3D B = (n_eo - (n_eo * e_b) * e_b) * (sin(KappaL) - KappaL) / Kappa_o;
    Vec3D C = (n_eo ^ e_b) * (1.0 - cos(KappaL)) / Kappa_o;
    DeltaR += (B + C);
  }

  Vec3D r_new = Lepton->m_r4.r + DeltaR;

  return r_new;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ICoord MagneticGrid::CheckCurrentCell(Vec3D r)
/*
  Calculates ICoord based on the current position in the grid. This is
  converts a VEC3D_T into a double, so it can then be converted to an
  integer to find the cell at a given particle position.

  \param  r   - position of the particle [cm]

  \return pos - ICoord, collection of 3 integers which
                denotes the grid cell.

*/
{

  double rx = Double(r.x);
  double ry = Double(r.y);
  double rz = Double(r.z);

  int n_x = (int)round(rx / m_cellsize);
  int n_y = (int)round(ry / m_cellsize);
  int n_z = (int)round(rz / m_cellsize);

  return ICoord(n_x, n_y, n_z);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CheckMagneticField
// \param: ICoord, an index for the map
inline Vec3D MagneticGrid::CheckMagneticField(const ICoord &c) {
  if (m_MagneticField.find(c) == m_MagneticField.end()) {
    // cell has not been used before so create the field
    // B bnew = generateRandomField();
    Vec3D bnew;
    bnew = bnew.UniformSphereDirection(m_rng);
    m_MagneticField[c] = bnew;
  }

  Vec3D b = m_MagneticField[c];

  return b;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline Vec3D MagneticGrid::CheckMagneticField_Lock(const ICoord &c)
/*
  Algorithm as follows:
    1) Check if ICoord defined
    2) If not, open file for reading, if ICoord in file, copy whole file
      into MagField map, return approp B value.
    3) If ICoord not in file, open file for reading, copy whole file into
      MagField Map, close for reading, lock immediately, open for writing,
      define B value at this ICoord, write map to file, lock/close.
*/

{

  if (m_MagneticField.find(c) == m_MagneticField.end()) {

    // char* filename = sfilename.c_str();

    ICoord c_ = c;

    // Open for reading, after obtaining the lock...
    std::ifstream infile;
    int *fd = new int;
    std::string rw_type = "read";
    int error_msg = LockAttempt(m_sfilename.c_str(), fd, rw_type);
    infile.open(m_sfilename.c_str());

    if (!infile) {
      std::cerr << "Unable to open: " << m_sfilename << std::endl;
      std::cerr << "ERROR CheckMagneticField_Lock 1." << std::endl;
      std::cerr << "error_msg = " << error_msg << std::endl;
      exit(EXIT_FAILURE); // call system to stop
    }

    int x = 0;
    int y = 0;
    int z = 0;
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;
    ICoord coord(x, y, z);
    Vec3D Bin(Bx, By, Bz);
    while (infile) {
      infile >> x >> y >> z;
      infile >> Bx >> By >> Bz;

      coord.ix = x;
      coord.iy = y;
      coord.iz = z;

      // If the grid cell has been defined already in the file
      if (c_ == coord) {

        ReadFieldMap(infile);
        Vec3D b = m_MagneticField[c];

        std::cerr << "ICoord found in MF file." << std::endl;

        infile.close();
        close(*fd);
        delete fd;

        return b;
      }
    }

    // Grid cell has not been defined in file.
    // Need: close file for reading, open for writing, AND lock! Then:
    //   1) read entire map from file into MagneticField[],
    //   2) define random field vector
    //   3) write map into text file

    infile.clear(); // Must clear ifstream before resetting file ptr!!
    infile.seekg(0, std::ios::beg);
    infile.open(m_sfilename.c_str());

    if (error_msg == -2) {
      std::cerr << "ERROR: Locking error 1." << std::endl;
      exit(EXIT_FAILURE);
    }

    ReadFieldMap(infile);

    infile.close();
    close(*fd);

    std::ofstream outfile;
    rw_type = "write";
    error_msg = LockAttempt(m_sfilename.c_str(), fd, rw_type);

    if (error_msg < 0) {
      std::cerr << "ERROR: Locking error 2." << std::endl;
      std::cerr << "error_msg = " << error_msg << std::endl;
      exit(EXIT_FAILURE);
    }

    outfile.open(m_sfilename.c_str());

    if (outfile.fail()) {
      std::cerr << "Unable to open: " << m_sfilename << std::endl;
      std::cerr << "ERROR CheckMagneticField_Lock 2." << std::endl;
      std::cerr << "error_msg = " << error_msg << std::endl;
      exit(EXIT_FAILURE);
    }

    Vec3D bnew;
    bnew = bnew.UniformSphereDirection(m_rng);
    m_MagneticField[c] = bnew;
    PrintFieldMapToFile(outfile);

    outfile.close();
    close(*fd);
    delete fd;

    return bnew;

  } else {
    Vec3D b = m_MagneticField[c];
    return b;
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  void MagneticGrid::ReadFieldMap(Field& MagneticField,
void MagneticGrid::ReadFieldMap(std::ifstream &infile)
/*

Reads data from filename textfile into m_MagneticField map.

*/
{

  infile.clear();
  infile.seekg(0, std::ios::beg);

  if (!infile) {
    std::cerr << "Unable to open infile." << std::endl;
    std::cerr << "ERROR ReadFieldMap." << std::endl;
    exit(EXIT_FAILURE); // call system to stop
  }

  int x = 0;
  int y = 0;
  int z = 0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  ICoord coord(x, y, z);
  Vec3D Bin(Bx, By, Bz);
  while (infile) {

    infile >> x >> y >> z;
    infile >> Bx >> By >> Bz;

    coord.ix = x;
    coord.iy = y;
    coord.iz = z;

    Bin.x = Bx;
    Bin.y = By;
    Bin.z = Bz;
    m_MagneticField[coord] = Bin;
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void MagneticGrid::PrintFieldMapToFile(std::ofstream &outfile)
/*
  Writes out m_MagneticField Map data into textfile filename.
*/
{

  if (!outfile) {
    std::cerr << "Unable to open outfile." << std::endl;
    std::cerr << "ERROR PrintFieldMapToFile." << std::endl;
    exit(EXIT_FAILURE);
  }

  for (Field::const_iterator it = m_MagneticField.begin();
       it != m_MagneticField.end(); ++it) {
    ICoord coord = it->first;
    Vec3D bvalue = it->second;
    outfile << coord << " ";
    outfile << bvalue.x << " " << bvalue.y << " " << bvalue.z << std::endl;
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int MagneticGrid::LockAttempt(const char *filename, int *fd, std::string rwtype)
/*
  file locking (so multiple users don't interfere with each other)

  \return: a file descriptor if successful, else negative number for
           error message.

  \note: locks are released when the files are closed in UnLock()
         function when file descriptor is closed.
*/
{

  const int WAIT_TIME = 10; // 10 musec
  struct flock lck;

  if (rwtype == "read") { // read_write = TRUE for read lock
    *fd = open(filename, O_RDONLY);
    lck.l_type = F_RDLCK;
    std::cerr << "Read lock attempted." << std::endl;
  } else if (rwtype == "write") {
    *fd = open(filename, O_WRONLY);
    lck.l_type = F_WRLCK;
    std::cerr << "Write lock attempted." << std::endl;
  } else {
    std::cerr << "Invalid Lock Type: " << rwtype << std::endl;
    exit(EXIT_FAILURE);
  }

  if (*fd < 0) return -1; // file not found, or could not be opened

  // set up the record locking structure, the address of which
  // is passed to the fcntl system call.
  lck.l_whence = SEEK_SET;
  lck.l_start = 0; // from beginning of file
  lck.l_len = 0;   // until end of the file

  // Attempt locking

  ///////////////////////////////////////////////////////////////
  // NOTE: errno codes are listed in /usr/include/errno.h
  //   Check them to see what the problem was if needed! (or online)
  ///////////////////////////////////////////////////////////////

  int TRYNO = 0;
  std::cerr << "Lock Attempt..." << std::endl;
  while (TRYNO >= 0) { // only a return can break out of the while loop

    if (fcntl(*fd, F_SETLK, &lck) == 0) {
      std::cerr << "TRYNO: " << TRYNO << std::endl;
      return *fd;
    } else if (errno == EAGAIN || errno == EACCES) { // file in use
      // std::cerr<<"File in use..."<<std::endl;
      usleep(WAIT_TIME);
    } else if (errno == EIO) {
      std::cerr << "IO Locking ERROR? errno: " << errno << std::endl;
      usleep(WAIT_TIME);
    } else {
      std::cerr << "Locking ERROR? errno: " << std::endl;
      std::cerr << "TRYNO: " << TRYNO << std::endl;
      usleep(WAIT_TIME);
      // std::cerr<<"errno: "<<errno<<std::endl;
      // return -2;
    }
    TRYNO = TRYNO + 1;
  }

  return -3; // nothing happened!
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} // namespace IGCascade
