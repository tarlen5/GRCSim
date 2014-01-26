/*!
-------------------------------------------------------------------------------
    \file   IGTrackElectrons.cpp

    Class that stores all variables related to intergalactic
    simulations. Also runs all qed processes.
  
    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
	       arlen@astro.ucla.edu                  \n

    \date      August 17, 2012                       \n
    
    \note 
-------------------------------------------------------------------------------
*/

#include "IGTrackElectrons.hpp"

using PhysConst::OMEGA_R; 
using PhysConst::OMEGA_M;
using PhysConst::OMEGA_L;
using PhysConst::OMEGA_0;

const bool SAVE_LOW_EGY = false;
const bool FORCE_SINGLE_GEN = false;
const bool TRACK_LEPTONS = false;

namespace IGCascade
{

  IGTrackElectrons::IGTrackElectrons(const string& s_egy, const string& s_Bmag, 
			     const string& s_ze, const string& s_cellsize, 
			     const string& s_iterations, 
			     const string& s_eblmodel,const string& s_file_num)
  {
    m_egy_cascade = s_egy.c_str();
    m_egy_cascade*=1.0e-3; // TeV
    m_bmag        = s_Bmag.c_str();
    m_ze          = s_ze.c_str();
    m_cellsize    = s_cellsize;
    double  d_iterations = 0.0;
    istringstream(s_iterations) >> d_iterations;
    m_iterations = (int) d_iterations;

    cout<<"Running sim with global variables:\n  SAVE_LOW_EGY = "<<SAVE_LOW_EGY<<"\n  FORCE_SINGLE_GEN = "<<FORCE_SINGLE_GEN<<"\n  TRACK_LEPTONS = "<<TRACK_LEPTONS<<endl;
    cout<<"Defining cascade file..."<<endl;
    m_cascade_file = DefineCascadeFile(s_eblmodel, s_egy, s_Bmag, s_ze, 
				       s_cellsize, s_file_num);

    string static_var_file = "StaticVariables.txt";
    DefineStaticVariables(static_var_file);
    DefineRp();

    m_rng = new RandomNumbers("random_numbers.seed");
    
    // NOTE: The first line of static_var_file MUST be filename for MFgrid
    string MFfilename = DefineMFfile(static_var_file,s_Bmag, s_cellsize, s_ze);

    m_BFieldGrid = 
      new MagneticGrid(m_rng, m_bmag, s_cellsize, MFfilename);
    m_pspace = new PairProduction(m_rng, m_ze);
    m_kspace = new KleinNishina(m_rng);

    if(SAVE_LOW_EGY) {
      m_low_egy_file = DefineLowEgyFile(s_eblmodel, s_egy, s_Bmag, s_ze,
					s_cellsize, s_file_num);
      ofstream low_egy_stream(m_low_egy_file.c_str());  // Overwrite existing.
      low_egy_stream.close(); 
    }

    if (TRACK_LEPTONS) {
      m_save_lepton_file = DefineTrackLeptonFile(s_eblmodel, s_egy, s_Bmag, 
						 s_ze,s_cellsize, s_file_num);
      ofstream ofile(m_save_lepton_file.c_str()); // Overwrite existing
      ofile.close();      
    }

    string ebl_model = s_eblmodel+".dat";
    cerr<<"File to be analyzed: "<<m_cascade_file<<endl;
    InitializeEBL(ebl_model);

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  string IGTrackElectrons::
  DefineCascadeFile(const string& s_eblmodel,const string& s_egy,
		    const string& s_Bmag, const string& s_ze, 
		    const string& s_cellsize, const string& s_file_num)
  {
    string filename = "TrackElectrons_"+s_eblmodel+"_"+s_egy+"GeV_B"+
      s_Bmag+"_z"+s_ze+"_L"+s_cellsize+"_"+s_file_num+".txt";
    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  string IGTrackElectrons::
  DefineMFfile(const string& static_var_file, const string& s_Bmag,
	       const string& s_cellsize, const string& s_ze)
  {
    string MFfilename;
    ifstream DefineFilename(static_var_file.c_str());
    DefineFilename >> MFfilename;
    DefineFilename.close();

    MFfilename+=s_Bmag+"_L"+s_cellsize+"_z"+s_ze+".txt";    
    cout<<"\nUsing Mag Grid file: "<<MFfilename<<endl;
    
    return MFfilename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  string IGTrackElectrons::
  DefineLowEgyFile(const string& s_eblmodel,const string& s_egy,
		   const string& s_Bmag, const string& s_ze, 
		   const string& s_cellsize, const string& s_file_num)
  {
    string filename = "LowEgy_"+s_eblmodel+"_"+s_egy+"GeV_B"+s_Bmag+"_z"+s_ze+
      "_L"+s_cellsize+"_"+s_file_num+".txt";
    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  string IGTrackElectrons::
  DefineTrackLeptonFile(const string& s_eblmodel, const string& s_egy, 
			const string& s_Bmag, const string& s_ze, 
			const string& s_cellsize, const string& s_file_num)
  {
    string filename = "TrackLeptons_"+s_eblmodel+"_"+s_egy+"GeV_B"+s_Bmag+"_z"+
      s_ze+"_L"+s_cellsize+"_"+s_file_num+".txt";
    
    return filename;
  }


  void IGTrackElectrons::DefineStaticVariables(const string& static_var_file)
  {
    
    ifstream ConstantList(static_var_file.c_str());
    if (ConstantList == NULL) {
      std::cerr << "could not open " << static_var_file << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Default values of static variables:
    int Max_Char = 500;
    char line[Max_Char];
    double temp_Tcmb = 0.0;
    double temp_E_min = 0.0;
    double temp_E_l_min = 0.0;
    bool temp_LOCK = 0;
    while ( !ConstantList.getline(line,Max_Char).eof() ) {      
      istringstream ss_line(line);
      string name;
      string svalue;
      ss_line >> name >> svalue;
      double value;
      istringstream istr(svalue);
      istr >> value;
      
      if (name[0] == '#') continue;          // skip over the commented line
      if (name == "Tcmb") temp_Tcmb = value;
      if (name == "E_min") temp_E_min = value;
      if (name == "E_l_min") temp_E_l_min = value;
      if (name == "LOCK") temp_LOCK = value;
    }
    const double Tcmb  = temp_Tcmb;            // [K]
    m_egy_gamma_min    = temp_E_min*PhysConst::TeV*0.99;  // TeV
    m_egy_lepton_min   = temp_E_l_min*PhysConst::TeV;     // TeV
    m_LOCK             = temp_LOCK;
    m_egy_cmb          = Tcmb*PhysConst::eV_K_B;
    
    // Testing:
    std::cout<<"Static Variables: "<<std::endl;
    std::cout<<"Tcmb = "<<Tcmb<<std::endl;
    std::cout<<"E_min = "<<m_egy_gamma_min<<std::endl;
    std::cout<<"E_l_min = "<<m_egy_lepton_min<<std::endl;
    std::cout<<"LOCK = "<<m_LOCK<<std::endl;
    std::cout.flush();
    ConstantList.close();
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    
        
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  void IGTrackElectrons::DefineRp(void)
  /*!
    Calculates magnitude from z=ze point to z_s = 0 surface, in order to come 
    up with a proper cutoff for integration for photons to the interaction
    redshift, in function PropagateDirectPhoton() and GetMinZ()

    NOTE: Instead of importing a file with the value, it doesn't take
    much time to merely calculate the value of Rp to very good
    precision. Values here computed with dz = 1.0e-5 were checked
    against the values in Rp_z_table.txt-which are the exact values of
    the integrals using the ARPREC code, and they all up to z=1.0 were
    exactly right to better than 1 part in 10^6!

    If we want to go past z=1, we should check with the exact value
    using the ARPREC code.

  */
  {
    
    /////// Now integrate it...
    double dz = 1.0e-5;
    double zfinal = Double(m_ze);
    double R_H = Double(PhysConst::CGS_HUBRAD);
    double OmegaR = Double(OMEGA_R);
    double OmegaM = Double(PhysConst::OMEGA_M);
    double OmegaL = Double(PhysConst::OMEGA_L);
    
    double zi = 0.0;
    double R_p = 0.0;
    unsigned num_steps = 0;
    while(zi < zfinal) {
      // zi is the low edge of the bin
      // z is the midpoint, where function is taken.
      
      double z = 0.0;
      if ((zi + dz) > zfinal) z = (zi + zfinal)/2.0;
      else  z = (zi + zi+dz)/2.0;
      
      double z1 = (1.0+z);
      double z2 = z1*z1;
      double z3 = z2*z1;
      double z4 = z3*z1;
      double Q = sqrt(OmegaR*z4 + OmegaM*z3 + OmegaL);
      
      if( (zi + dz) > zfinal) {
	R_p += (zfinal-zi)/Q;
	zi = zfinal;
      }
      else {
	R_p += dz/Q;
	zi+=dz;
      }
      num_steps++;
    }
    //////////////////////////////////////
    
    m_R_0 = R_p*R_H;

    // cout<<"R_p from integral: "<<R_p<<"\n  in "<<num_steps<<" steps."<<endl;
    // cout<<"percent error: "<<fabs(R_p - m_R_0)/m_R_0*100.0<<endl;
    // char getline;
    // cin>>getline;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
  

  void IGTrackElectrons::InitializeEBL(const string& ebl_model_file)
  {

    // The string "table" in the filename is what will distinguish the
    // table file from the normal DIRBR file.
    if( ebl_model_file.find("table") == string::npos ) {

      m_ebl = new DIRBR();
      
      //////////////////////////////////////////////////////////////
      //-------------------INITIALIZE EBL MODEL---------------------
      //////////////////////////////////////////////////////////////
      std::ifstream model;
      model.open(ebl_model_file.c_str());
      
      if (model == NULL) {
	std::cerr << "ERROR: could not open EBL model file: " << ebl_model_file 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      
      vector<double> lambda_vec;
      vector<double> nuFnu_vec;
      // Set up DIRBR curve...
      while(model) {
	double lambda_i;
	double nuFnu_i;
	model >> lambda_i >> nuFnu_i;
	if(model) {
	  lambda_vec.push_back(lambda_i);
	  nuFnu_vec.push_back(nuFnu_i);
	}
      }
      model.close();

      cout<<"initializing DIRBR..."<<endl;

      // Initialize DIRBR:
      DIRBR_ERR ebl_err=m_ebl->SetDIRBR(lambda_vec.size(), &lambda_vec[0], 
				       &nuFnu_vec[0]);
      std::cout << "SetDIRBR returned: " << ebl_err << std::endl;
      ebl_err=m_ebl->TestDIRBRlimits();
      std::cout << "TestDIRBRlimits returned: " << ebl_err << std::endl;
    
    } else {
      
      // Code for declaring the 2D table version of DIRBR EBL Model.
      
    }
    
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  void IGTrackElectrons::RunCascade(void)
  {

    // -----------Print time at beginning-----------
    time_t curr=time(0);
    std::cerr<<"Start Time: "<<ctime(&curr)<<std::endl<<std::endl;
    //----------------------------------------------
    
    ofstream PhotonList(m_cascade_file.c_str()); // Overwrite existing
    PhotonList.close();
    
    //  int max_iterations = 10;
    for (int i=1; i<=m_iterations; i++) {
      
      // This is printed out every time a new source photon is defined
      PhotonList.open(m_cascade_file.c_str(), std::ios::app);
      PhotonList<<endl<<m_ze<<" "<<m_BFieldGrid->m_bmag<<" "<<
	m_BFieldGrid->m_cellsize/PhysConst::MPC_TO_CM<<" "<<m_egy_cascade<<endl;
      PhotonList.close();

      RelParticle GammaPhoton;
      CreateDirectGamma(GammaPhoton);
      
#ifdef DEBUG  
      cout<<"Photon iteration: "<<i<<endl;
      cout<<endl<<m_ze<<" "<<m_BFieldGrid->m_bmag<<" [gauss] "<<
	m_BFieldGrid->m_cellsize/PhysConst::MPC_TO_CM<<" [Mpc] "<<m_egy_cascade
	  <<" [TeV]"<<endl;
      cout<<"\nGamma created with parameters [TeV]: \n  egy: "<<
	GammaPhoton.m_p4.r0<<"\n  3 momentum: "<<GammaPhoton.m_p4.r<<endl;
#endif // DEBUG

      // Cascade initiated from Direct photon:
      RelParticle* Electron = new RelParticle;
      Electron->Zero();
      RelParticle* Positron = new RelParticle;
      Positron->Zero();
      bool pair_prod = PropagateDirectPhoton(GammaPhoton,Electron,Positron);
      
      stack<RelParticle*> lepton_stack;
      unsigned stack_max = 0;
      
      PhotonList.open(m_cascade_file.c_str(), std::ios::app);
      
      unsigned ilepton = 0;
      if(!pair_prod) SaveDirectPhoton(GammaPhoton, PhotonList);
      else {  // Cascaded Leptons initiated...
	
	GammaPhoton.Zero();

	if(Positron->m_p4.r0 > m_egy_lepton_min) {
	  lepton_stack.push(Positron);
	  ilepton++;
	  Positron->m_tag = ilepton;
	  if (TRACK_LEPTONS) SaveLepton(Positron);
	}
	else {
	  if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Positron);
	  delete Positron;
	}
	if(Electron->m_p4.r0 > m_egy_lepton_min) {
	  lepton_stack.push(Electron);
	  ilepton++;
	  Electron->m_tag = ilepton;
	  if (TRACK_LEPTONS) SaveLepton(Electron);
	}
	else {
	  if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Electron);
	  delete Electron;
	}

	// exit condition is when there are no more leptons left for
	// IC Scattering
	stack_max = lepton_stack.size();
	while( !lepton_stack.empty() ) {
	  
	  RelParticle* Lepton = lepton_stack.top();
	  lepton_stack.pop();
	  PropagateLepton(Lepton, GammaPhoton);
	  
	  if(Lepton->m_p4.r0 > m_egy_lepton_min && 
	     Lepton->m_z_s > m_BFieldGrid->m_DE) {
	    lepton_stack.push(Lepton);
	    if (TRACK_LEPTONS) SaveLepton(Lepton);
	  } else {
	    if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Lepton);
	    cerr<<"lepton z_s = "<<Lepton->m_z_s<<" egy: "<<Lepton->m_p4.r0
		<<endl;
	    delete Lepton;
	    cerr<<"lepton deleted from stack"<<endl<<endl;
	  }
	  
	  //////////////////////////////////////////////////////////////////
	  // IF GAMMA ENERGY TOO LOW OF IC PHOTON BEFORE PROPAGATING IT...
	  // Don't bother propagating to z_s = 0 surface, just reject it now.
	  //////////////////////////////////////////////////////////////////
	  if ( GammaPhoton.m_p4.r0 < m_egy_gamma_min ) {
#ifdef DEBUG
	    cout<<"  IC Gamma egy too low, egy: "<<GammaPhoton.m_p4.r0<<endl;
#endif
	    if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
	    GammaPhoton.Zero();
	  } else {
	    // Propagate Secondary Gamma...
	    RelParticle* Electron = new RelParticle;
	    Electron->Zero();
	    RelParticle* Positron = new RelParticle;
	    Positron->Zero();
	    
	    bool pair_prod = 
	      PropagateSecondaryPhoton(GammaPhoton,Electron,Positron);
	    
	    if (!pair_prod) {
	      if(GammaPhoton.m_z_s > 1.0e-15) { // Should never happen
		cerr<<"\nWARNING: Caught z_s=0 photon with z_s> 1E-15!"<<endl;
		cerr<<"(Secondary photon...)"<<endl<<endl;
	      }
	      
	      if (GammaPhoton.m_p4.r0 >= m_egy_gamma_min)
		m_pspace->WritePhotonToFile(PhotonList, GammaPhoton.m_p4, 
				       GammaPhoton.m_r4,GammaPhoton.m_weight);
	      else {
		if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
	      }
	    }
	    else { // Add electrons to stack?
	      if(Positron->m_p4.r0 >m_egy_lepton_min) {
		lepton_stack.push(Positron);
		ilepton++;
		Positron->m_tag = ilepton;
		if (TRACK_LEPTONS) SaveLepton(Positron);
	      }
	      else {
		if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Positron);
		delete Positron;
	      }
	      if(Electron->m_p4.r0>m_egy_lepton_min) {
		lepton_stack.push(Electron);
		ilepton++;
		Electron->m_tag = ilepton;
		if (TRACK_LEPTONS) SaveLepton(Electron);
	      }
	      else {
		if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Electron);
		delete Electron;
	      }
	    }
	    if(lepton_stack.size() > stack_max) stack_max = lepton_stack.size();
	    GammaPhoton.Zero();
	  }
	  
	} // End while loop over leptons
      } // End if cascade initialized
      PhotonList.close();
      cerr<<"\nMax number of leptons in stack = "<<stack_max<<endl<<endl;

    } // End loop over iterations
      
    //----------Print End Time----------
    curr=time(0);
    std::cerr<<"End Time: "<<ctime(&curr)<<std::endl;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  void IGTrackElectrons::CreateDirectGamma(RelParticle& GammaPhoton)
  { 

    GammaPhoton.Zero();
    
    //---------Initialize Gamma Photon along z-axis---------//
    VEC3D_T E_e = m_egy_cascade*PhysConst::TeV*(1.0+m_ze);   // energy at EMISSION
    Vec4D R4_e(0.,0.,0.,0.);
    Vec3D n1(0.,0.,1.);
    GammaPhoton.m_p4.r0 = E_e;
    GammaPhoton.m_m0 = 0.0;
    GammaPhoton.m_r4 = R4_e;
    GammaPhoton.m_p4.r = n1*E_e;
    GammaPhoton.m_z = m_ze;
    GammaPhoton.m_z_s = m_ze;
    GammaPhoton.m_egy_int = E_e;
    GammaPhoton.m_z_s_int = "0.0";
    GammaPhoton.m_z_int = "0.0";
    GammaPhoton.m_rnext = "0.0";
    GammaPhoton.m_q = 0;
    GammaPhoton.m_weight = 1.0;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  bool IGTrackElectrons::
    PropagateDirectPhoton(RelParticle& GammaPhoton, RelParticle* Electron, 
			  RelParticle* Positron)
  /*!  
    At the end of this function, the Electron and Positron are
    created at the interaction redshift, IF a pair production event
    takes place.
  */
  {
    
    double z_min = GetMinZ(GammaPhoton.m_p4, GammaPhoton.m_r4, GammaPhoton.m_z);
#ifdef DEBUG 
    cout<<"\nPropagating DIRECT photon through EBL...\n";
    cout<<"z_min = "<<z_min<<endl;
#endif
    double tot_lambda_int = 0.0;
    VEC3D_T z_int = "0.0";
    VEC3D_T z_emit = GammaPhoton.m_z;
    bool pair_prod =  
      m_pspace->PropagatePhotonEBL(m_ebl,Double(GammaPhoton.m_p4.r0),
				   Double(GammaPhoton.m_z),z_min,z_int,
				   tot_lambda_int);
    
    if(!pair_prod) {
      cout<<"NO PAIR PRODUCTION...\n";
      // Propagate to z_s = 0 surface...
      //GammaPhoton.m_p4.r0 = 0.0;
      Vec3D Rfinal(0.0,0.0,m_R_0);
      GammaPhoton.m_r4.r  = Rfinal;
      GammaPhoton.m_p4.r0 *= 1.0/(1.0+GammaPhoton.m_z);
      GammaPhoton.m_p4.r  *= GammaPhoton.m_p4.r0/GammaPhoton.m_p4.r.Norm();
      GammaPhoton.m_z     = 0.0;
      GammaPhoton.m_z_s   = 0.0;
      
      return false;
    }
    // z_min is not reached...However, from propagating photon
    // to z_final, it is possible that z_s=0 surface could be
    // crossed, since z_min does not correspond exactly to the
    // z_s=0 surface for computation time issues.
    else { 
      VEC3D_T delta_z = (z_emit - z_int);
      bool pair_prod = 
	m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				    GammaPhoton.m_z,GammaPhoton.m_z_s, 
				    delta_z);
      if( !pair_prod ) return false;

      // Otherwise, pair production has occured and proceed as usual
      // defining EBL photon interacted with, Leptons, etc.
      RelParticle EBLPhoton;
      EBLPhoton.m_p4.r0 = m_pspace->GetEBLPhotonEgy(m_ebl,tot_lambda_int,
					   Double(z_int),Double(z_emit),
					   Double(GammaPhoton.m_p4.r0));
      m_pspace->UpdateEBLPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				GammaPhoton.m_z,GammaPhoton.m_z_s,
				EBLPhoton.m_p4,EBLPhoton.m_r4,EBLPhoton.m_z);
      
      ////////////////////////////////////////////////////
      // CREATE AND DEFINE LEPTONS PRODUCED BY PAIRPROD //
      ////////////////////////////////////////////////////
      DefineLeptons(GammaPhoton,EBLPhoton,Electron,Positron);     

      return true;
    }
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGTrackElectrons::
  DefineLeptons(RelParticle& GammaPhoton, RelParticle& EBLPhoton, RelParticle*
		Electron, RelParticle* Positron)
  /*!  
    Uses GammaPhoton and EBLPhoton to define Electron and Positron
    in the lab frame.
    
   */
  {
    
    ////////////////////////////////////////////////////
    // CREATE AND DEFINE LEPTONS PRODUCED BY PAIRPROD //
    ////////////////////////////////////////////////////
    std::cerr<<"Pair Production Occurs..."<<std::endl;
    int charge_elec = -1;
    Vec4D p4(PhysConst::SI_MELEC,0.0,0.0,0.0);
    Electron->m_r4      = GammaPhoton.m_r4;
    Electron->m_p4      = p4;
    Electron->m_m0      = PhysConst::SI_MELEC;
    Electron->m_z       = GammaPhoton.m_z;
    Electron->m_z_s     = GammaPhoton.m_z_s;
    Electron->m_z_s_int = GammaPhoton.m_z_s;
    Electron->m_z_int   = GammaPhoton.m_z;
    Electron->m_rnext   = GammaPhoton.m_rnext;
    Electron->m_q       = charge_elec;
    Electron->m_weight  = GammaPhoton.m_weight;
    std::cerr<<endl<<"new RelParticle: Electron..."<<endl;
    
    Positron->m_r4 =  Electron->m_r4;
    Positron->m_p4 =  Electron->m_p4;
    Positron->m_m0 =  Electron->m_m0;
    Positron->m_z  =  Electron->m_z;
    Positron->m_z_s = Electron->m_z_s;
    Positron->m_z_s_int = GammaPhoton.m_z_s;
    Positron->m_z_int = GammaPhoton.m_z;
    Positron->m_rnext = GammaPhoton.m_rnext;
    Positron->m_q  = -Electron->m_q;
    Positron->m_weight = Electron->m_weight;
    std::cerr<<endl<<"new RelParticle: Positron..."<<endl;
    /////////////////////////////////////////
    
    bool RelKinematics = 
      m_pspace->RelativisticKinematics(GammaPhoton.m_p4,EBLPhoton.m_p4,
				       Electron->m_p4, Positron->m_p4);
    
    if (!RelKinematics)
      std::cerr<<std::endl<<"ERROR: RelativisticKinematics failed...\n\n";
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  bool IGTrackElectrons::
  PropagateSecondaryPhoton(RelParticle& GammaPhoton, RelParticle* Electron, 
			   RelParticle* Positron)
  {

    //#ifdef DEBUG 
    cout<<"Propagating SECONDARY photon through EBL...";
    //#endif
    double z_min = GetMinZ(GammaPhoton.m_p4, GammaPhoton.m_r4, GammaPhoton.m_z);
    VEC3D_T z_emit = GammaPhoton.m_z;

    if (FORCE_SINGLE_GEN) {
      VEC3D_T delta_z = (z_emit - z_min);
      bool pair_prod = m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,
						   GammaPhoton.m_r4,
				  GammaPhoton.m_z,GammaPhoton.m_z_s, 
				  delta_z);
      return false;
    }
    
    double tot_lambda_int = 0.0;
    VEC3D_T z_int = "0.0";
    bool pair_prod =  
      m_pspace->PropagatePhotonEBL(m_ebl,Double(GammaPhoton.m_p4.r0),
				   Double(z_emit),z_min,z_int,
				   tot_lambda_int);
    // At this point, either z_min is reached so that no pair
    // production will occur, or z_int is reached at which pair
    // production happens. Either way, propagate gamma photon over
    // delta_z, and see if it reaches the z_s = 0 surface.

    //#ifdef DEBUG
    //cout<<"Secondary Photon propagation complete...Parameters:\n";
    //cout<<"egy: "<<GammaPhoton.m_p4.r0<<" zs: "<<GammaPhoton.m_z_s<<endl;
    //#endif        

    VEC3D_T delta_z = (z_emit - z_int);
    pair_prod = 
      m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				  GammaPhoton.m_z,GammaPhoton.m_z_s, 
				  delta_z);
    if( !pair_prod ) return false;
    // Otherwise, pair production has occured and proceed as usual
    // defining EBL photon interacted with, Leptons, etc.
    RelParticle EBLPhoton;
    EBLPhoton.m_p4.r0 = m_pspace->GetEBLPhotonEgy(m_ebl,tot_lambda_int,
						   Double(z_int),Double(z_emit),
						   Double(GammaPhoton.m_p4.r0));
    m_pspace->UpdateEBLPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
			      GammaPhoton.m_z,GammaPhoton.m_z_s,
			      EBLPhoton.m_p4,EBLPhoton.m_r4,EBLPhoton.m_z);
    
    ////////////////////////////////////////////////////
    // CREATE AND DEFINE LEPTONS PRODUCED BY PAIRPROD //
    ////////////////////////////////////////////////////
    DefineLeptons(GammaPhoton,EBLPhoton,Electron,Positron);     
    
    return true;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
  double IGTrackElectrons::GetMinZ(const Vec4D& gam_ph_p4, const Vec4D& gam_ph_r4,
			       const VEC3D_T& gam_ph_z)
  /*!  
    Quick routine to compute z_min, the change in redshift (not
    change in z_s, change in z) to get to the z_s = 0 surface.

    WARNING: This OVERESTIMATES |z_min|, so that if the true z_min
    were e.g. -0.001 it would estimate at something like -0.00115, so
    that in the numerical integration for tau, there will be a buffer
    to go past this value a little and guarantee that it will reach
    the z_s = 0 surface when required .
    
    \param gam_ph_z - photon's redshift at beginning of propagation.

  */
  {

    ///////////////////////////////////////////////
    // Calculate z_min:
    //////////////////////////////////////////////
    VEC3D_T R_test = m_R_0;
    VEC3D_T R1 = gam_ph_r4.r.Norm();
    Vec3D e2 =  gam_ph_p4.r/gam_ph_p4.r.Norm();
    Vec3D e1(e2);
    if ( R1 > 0.0 ) e1 = gam_ph_r4.r/gam_ph_r4.r.Norm();    
    VEC3D_T costheta = e1*e2;
    VEC3D_T sintheta_sq = 1.0 - costheta*costheta;
    // Dist photon must travel for z_s = 0 surface:
    VEC3D_T R2 = sqrt( R_test*R_test - R1*R1*sintheta_sq) - R1*costheta;
    
    // Checking:
    Vec3D r_test = R1*e1 + R2*e2;
    R_test = r_test.Norm();
    
    // Overestimate z_min:
    VEC3D_T z1 = (1.0+gam_ph_z);
    VEC3D_T Q_hp = sqrt(OMEGA_R*z1*z1*z1*z1 + OMEGA_M*z1*z1*z1 + OMEGA_L +
		     (1.0-OMEGA_0)*z1*z1);
    // gives us z_min; must be smaller than z1 - DeltaZ.
    VEC3D_T DeltaZ = R2/PhysConst::CGS_HUBRAD*Q_hp;
    VEC3D_T z_min = gam_ph_z - DeltaZ;
    VEC3D_T num_int_steps = "200.0";
    VEC3D_T z_step = (gam_ph_z - z_min)/num_int_steps;
    
    R_test = "0.0";
    VEC3D_T zi = gam_ph_z - z_step/2.0;  // Midway between z[0] & z[1].
    int steps = 0;
    // integrate using simple simpson's midpoint method
    //std::cout<<"m_R_0: "<<m_R_0<<" R2: "<<R2<<std::endl;
    while(R_test < R2) {  // integrate until you reach the distance of R2.
      zi -= z_step;
      z1 = (1.0 + zi);
      Q_hp = sqrt(OMEGA_R*z1*z1*z1*z1 + OMEGA_M*z1*z1*z1 + OMEGA_L +
		  (1.0-OMEGA_0)*z1*z1);
      
      R_test += PhysConst::CGS_HUBRAD/Q_hp*z_step;
      steps++;
    }

    return Double(zi);
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGTrackElectrons::PropagateLepton(RelParticle* Lepton, 
				     RelParticle& GammaPhoton)
  {

#ifdef DEBUG
    cout<<"Propagating Lepton, at z = "<<Lepton->m_z_s<<" egy: "<<
      Lepton->m_p4.r0<<endl;
#endif

    VEC3D_T Ecmb_e = m_egy_cmb*(1.0+Lepton->m_z);
    VEC3D_T PL=m_kspace->PropagationLengthCMBandIR(m_ebl, Ecmb_e,Lepton->m_p4,
						GammaPhoton.m_p4);   // [cm]
    
    GammaPhoton.m_weight = Lepton->m_weight;
    
    VEC3D_T delta_z = m_BFieldGrid->Delta_z(PL, Lepton);
    
    // Lepton direction at zi, before scattering kinematics computed
    Vec3D n_e = Lepton->m_p4.r/Lepton->m_p4.r.Norm();
    m_kspace->RelativisticKinematics(Lepton->m_p4,GammaPhoton.m_p4);
    
    m_BFieldGrid->PropagateBFieldRedshift(GammaPhoton,Lepton,n_e,
					 PL,delta_z, m_LOCK);
    
    // NOTE: added to keep track of gamma_photon's emitted z_s
    GammaPhoton.m_z_s_int = Lepton->m_z_s;
    GammaPhoton.m_z_int = Lepton->m_z;
    GammaPhoton.m_egy_int = GammaPhoton.m_p4.r0;
    GammaPhoton.m_rnext = Lepton->m_rnext;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGTrackElectrons::
  SaveDirectPhoton(RelParticle& GammaPhoton, ofstream& PhotonList)
  {
        
    if(GammaPhoton.m_z_s > 1.0e-15) { // Should never happen
      cerr<<"\nWARNING: Caught z_s=0 photon with z_s > 1E-15!\n"<<endl;
      cerr<<"z_s: "<<GammaPhoton.m_z_s<<" egy: "<<GammaPhoton.m_p4.r0<<endl;
    }
    
    if (GammaPhoton.m_p4.r0 >= m_egy_gamma_min) {
      m_pspace->WritePhotonToFile(PhotonList, GammaPhoton.m_p4, 
				  GammaPhoton.m_r4, GammaPhoton.m_weight);
    }
    else if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
    

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGTrackElectrons::SaveToLowEnergyFile(RelParticle& Particle)
  {
    
    ofstream low_egy_stream(m_low_egy_file.c_str(), std::ios::app);
    low_egy_stream<<Particle.m_p4.r0<<" "<<Particle.m_z<<" "<<Particle.m_q<<
      endl;
    low_egy_stream.close();
    // Nothing yet...    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  void IGTrackElectrons::SaveLepton(RelParticle* Lepton)
  {
    
    // Prints out the line:
    //   number_tag z E px py pz time x y z
    //
    ofstream ofile(m_save_lepton_file.c_str(), std::ios::app);
    ofile<<Lepton->m_tag<<" "<<Lepton->m_z<<" "<<Lepton->m_p4.r0<<" "<<Lepton->m_p4<<" "<<Lepton->m_r4.r0<<" "<<Lepton->m_r4.r<<endl;
    ofile.close();

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
}

