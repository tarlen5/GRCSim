/*!
-------------------------------------------------------------------------------
    \file   IGCascadeSim.cpp

    Class than stores all variables related to intergalactic
    simulations. Also runs all qed processes.
  
    \author    Timothy C. Arlen                      \n
               Department of Physics and Astronomy   \n
               UCLA                                  \n
	       arlen@astro.ucla.edu                  \n

    \date      May 20, 2012                          \n
    
    \revision  v1.1  April 19, 2013
               Fairly major revision, by replacing all output to text
               files with all output to ROOT files, using mostly
               TTree.  This also requires rewrites of the intermediate
               scripts to convert the raw output from teh simulatios
               into the workable rootfiles per model of EBL,z,B,Lcoh.
-------------------------------------------------------------------------------
*/

#include "IGCascadeSim.hpp"

using PhysConst::OMEGA_R; 
using PhysConst::OMEGA_M;
using PhysConst::OMEGA_L;
using PhysConst::OMEGA_0;


namespace IGCascade
{

  IGCascadeSim::
  IGCascadeSim(const string& egy, const string& redshift, 
	       const string& mag_field, const string& coh_len, 
	       const string& file_count, AnyOption* opt)
  {
    
    string eblmodel, mf_dir, opt_depth_dir, output_dir;
    ProcessOptions(opt,eblmodel,mf_dir,opt_depth_dir,output_dir);

    string ebl_model_file = eblmodel+".dat";
    InitializeEBL(ebl_model_file);
    m_cascade_file = DefineCascadeFile(eblmodel, egy, mag_field, redshift, 
				       coh_len, file_count, output_dir);
    string MFfilename = DefineMFfile(mf_dir, mag_field, coh_len, redshift);
    string optDepthFile = DefineOptDepthTable(opt_depth_dir,eblmodel,redshift);
    
    if (m_trk_leptons_bool)
      m_save_lepton_file = DefineTrackLeptonFile(eblmodel, egy, mag_field,
				    redshift,coh_len, file_count, output_dir);
    
    if (m_trk_delay_bool)
      m_track_time_delay_file = 
	DefineTrackTimeDelayFile(eblmodel,egy,mag_field,redshift,coh_len,
				 file_count);
    
    // Process Inputs to class:
    m_egy_cascade = egy.c_str();
    m_egy_cascade*=1.0e-3;       // TeV
    m_ze          = redshift.c_str();
    m_bmag        = mag_field.c_str();
    m_cellsize    = coh_len.c_str();    
    
    m_DE = "1.0e-25";
    VEC3D_T Tcmb = 2.73;
    m_egy_cmb = Tcmb*PhysConst::eV_K_B;

    DefineRp();
    
    int file_num = atoi(file_count.c_str());
    m_rng = new TRandom3(0);
    m_rng->SetSeed(file_num);
    m_BFieldGrid =
      new MagneticGrid(m_rng, m_bmag, m_cellsize, MFfilename);
    m_pspace = new PairProduction(m_rng, m_ze);
    m_kspace = new KleinNishina(m_rng);

    cout<<endl;
    cout<<"Using cascade file:   "+m_cascade_file<<endl;
    cout<<"Using Mag Grid file:  "<<MFfilename<<endl;
    cout<<"Using opt depth file: "<<optDepthFile<<endl;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGCascadeSim::ProcessOptions(AnyOption* opt, string& eblmodel,
				    string& mf_dir, string& opt_depth_dir,
				    string& output_dir)
  {
    // Process Options:
    eblmodel      = "EBLModel4msld";
    mf_dir        = "MagneticFieldFiles/";
    opt_depth_dir = "OptDepthFiles/";
    output_dir    = "SimOutputFiles/";
    m_egy_gamma_min = 0.1;
    m_egy_lepton_min = 75.0;
    if(opt->getValue("eblmodel") != NULL) eblmodel = opt->getValue("eblmodel");
    if(opt->getValue("mf_dir") != NULL) mf_dir = opt->getValue("mf_dir");
    if(opt->getValue("opt_depth_dir") != NULL) 
      opt_depth_dir = opt->getValue("opt_depth_dir");
    if(opt->getValue("output_dir") != NULL) 
      output_dir = opt->getValue("output_dir");
    if(opt->getValue("gam_egy_min") != NULL)
      m_egy_gamma_min = opt->getValue("gam_egy_min");
    if(opt->getValue("lep_egy_min") != NULL)
      m_egy_lepton_min = opt->getValue("lep_egy_min");
    
    // Convert to eV:
    m_egy_gamma_min *= 1.0e9*0.99;
    m_egy_lepton_min *= 1.0e9;

    // Flags:
    m_LOCK = true;
    m_single_gen_bool = false;
    m_trk_delay_bool = false;
    m_trk_leptons_bool = false;
    if(opt->getFlag("mf_no_lock") != NULL)
      m_LOCK = false;
    if(opt->getFlag("single_gen") != NULL)
      m_single_gen_bool = true;
    if(opt->getFlag("trk_delay") != NULL)
      m_trk_delay_bool = true;
    if(opt->getFlag("trk_leptons") != NULL)
      m_trk_leptons_bool = true;

    cout<<"\n>> Options processed"<<endl;
    cout<<"   --eblmodel:      "<<eblmodel<<endl;
    cout<<"   --mf_dir:        "<<mf_dir<<endl;
    cout<<"   --opt_depth_dir: "<<opt_depth_dir<<endl;
    cout<<"   --output_dir:    "<<output_dir<<endl;
    cout<<"   --gam_egy_min:   "<<m_egy_gamma_min<<endl;
    cout<<"   --lep_egy_min:   "<<m_egy_lepton_min<<endl;
    cout<<"   --mf_no_lock:    "<<(!m_LOCK)<<endl;
    cout<<"   --single_gen:    "<<m_single_gen_bool<<endl;
    cout<<"   --trk_delay:     "<<m_trk_delay_bool<<endl;
    cout<<"   --trk_leptons:   "<<m_trk_leptons_bool<<endl;
    cout<<"-------------------------------------"<<endl<<endl;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  string IGCascadeSim::
  DefineCascadeFile(const string& s_eblmodel,const string& s_egy,
		    const string& s_Bmag, const string& s_ze, 
		    const string& s_cellsize, const string& s_file_num,
		    const string& output_dir)
  {
    string filename = output_dir+s_eblmodel+"_"+s_egy+"GeV_z"+s_ze+"_B"+s_Bmag+"_L"+
      s_cellsize+"_"+s_file_num+".root";
    
    m_secPhotonTree = 
      new TTree("Secondary","egyPrim,egySec,theta,phi,time,thetap,xi,weight");
    m_secPhotonTree->Branch("egyPrim",&m_egyPrim,"egyPrim/D");
    m_secPhotonTree->Branch("egySec",&m_egySec,"egySec/D");
    m_secPhotonTree->Branch("theta",&m_theta,"theta/D");
    m_secPhotonTree->Branch("phi",&m_phi,"phi/D");
    m_secPhotonTree->Branch("time",&m_time,"time/D");
    m_secPhotonTree->Branch("thetap",&m_thetap,"thetap/D");
    m_secPhotonTree->Branch("xi",&m_xi,"xi/D");
    m_secPhotonTree->Branch("weight",&m_weight,"weight/D");
    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  string IGCascadeSim::
  DefineMFfile(const string& mf_dir, const string& s_Bmag,
	       const string& s_cellsize, const string& s_ze)
  {

    string MFfilename=mf_dir+"MagneticGrid_B"+s_Bmag+"_L"+s_cellsize+"_z"+s_ze+".txt";    
    return MFfilename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  string IGCascadeSim::
  DefineLowEgyFile(const string& s_eblmodel,const string& s_egy,
		   const string& s_Bmag, const string& s_ze, 
		   const string& s_cellsize, const string& s_file_num)
  {
    string filename = "LowEgy_"+s_eblmodel+"_"+s_egy+"GeV_z"+s_ze+"_B"+s_Bmag+
      "_L"+s_cellsize+"_"+s_file_num+".txt";    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  string IGCascadeSim::
  DefineTrackLeptonFile(const string& s_eblmodel, const string& s_egy, 
			const string& s_Bmag, const string& s_ze, 
			const string& s_cellsize, const string& s_file_num,
			const string& output_dir)
  {
    string filename = output_dir+"TrackLeptons_"+s_eblmodel+"_"+s_egy+"GeV_z"+
      s_ze+"_B"+s_Bmag+"_L"+s_cellsize+"_"+s_file_num+".root";

    // Overwrite existing
    TFile* rootfile = new TFile(filename.c_str(),"RECREATE");
    rootfile->Close();
    ofstream ofile(m_track_time_delay_file.c_str()); // Overwrite existing
    ofile.close();
    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  string IGCascadeSim::
  DefineTrackTimeDelayFile(const string& s_eblmodel, const string& s_egy, 
			const string& s_Bmag, const string& s_ze, 
			const string& s_cellsize, const string& s_file_num)
  {
    string filename = "TrackTimeDelay_"+s_eblmodel+"_"+s_egy+"GeV_B"+s_Bmag+"_z"+
      s_ze+"_L"+s_cellsize+"_"+s_file_num+".txt";
    
    return filename;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  void IGCascadeSim::DefineRp(void)
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
 
  

  void IGCascadeSim::InitializeEBL(const string& ebl_model_file)
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

  string IGCascadeSim::
  DefineOptDepthTable(const string& opt_depth_dir, const string& s_eblmodel, 
		      const string& s_ze)
  {
    string optDepthFile = opt_depth_dir;
    optDepthFile+="optDepth_"+s_eblmodel+"_z"+s_ze+"_0.01TeV_1TeV.txt";
    
    m_optDepthTable = new Table2D(optDepthFile);
    m_tauCutoff = 0.003;  // This is the approx tau for a 56 GeV
			  // photon at z = 0.1
    
    return optDepthFile;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  void IGCascadeSim::RunCascade(const int numIterations)
  {

    // -----------Print time at beginning-----------
    time_t curr=time(0);
    std::cerr<<"\nStart Time: "<<ctime(&curr)<<std::endl<<std::endl;
    //----------------------------------------------

    // Overwrite Existing:
    TFile* rootfilePhotonList = new TFile(m_cascade_file.c_str(),"RECREATE");
    rootfilePhotonList->Close();

    m_globalLeptonNum = 0;
    for (int i=1; i<=numIterations; i++) {
      cout<<"-----------------------------------------------------------------";
      cout<<"\nNew Photon Defined, iteration = "<<i<<endl;
      cout<<"  egy: "<<m_egy_cascade<<" TeV"<<endl;

      RunSinglePhotonCascade();

    } // End loop over iterations

    // Saving Tree to rootfile:
    rootfilePhotonList = new TFile(m_cascade_file.c_str(),"UPDATE");
    m_secPhotonTree->Write();
    rootfilePhotonList->Close();
    delete rootfilePhotonList;
      
    //----------Print End Time----------
    double seconds = difftime(time(0),curr);
    std::cout<<"\nSimulation took: "<< (int(seconds)/3600)<< " hr, "<<
      (int(seconds)%3600)/60<<" min, "<< int(seconds)%60<<" sec."<<endl;
    cout<<"seconds: "<<seconds<<endl;
    curr=time(0);
    std::cerr<<"End Time: "<<ctime(&curr)<<std::endl;
      
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGCascadeSim::RunSinglePhotonCascade(void)
  {

    RelParticle GammaPhoton;
    CreateDirectGamma(GammaPhoton);
    
    stack<RelParticle*> lepton_stack;
    unsigned stack_max = 0;    
    cout<<"Propagating Direct Photon..."<<endl;
    //bool pair_prod = PropagateDirectPhoton(GammaPhoton,Electron,Positron);    
    PropagateDirectPhoton(GammaPhoton,lepton_stack);
    
    stack_max = lepton_stack.size();
    while( !lepton_stack.empty() ) {
      
      //cout<<"Propagating lepton..."<<endl;
      PropagateLepton(lepton_stack, GammaPhoton);
      //cout<<"Finished..."<<endl;
      
      //cout<<"  IC scattered photon created with egy: "<<GammaPhoton.m_p4.r0<<
      //" z_s: "<<GammaPhoton.m_z_s<<endl;
      
      //////////////////////////////////////////////////////////////////
      // IF GAMMA ENERGY TOO LOW OF IC PHOTON BEFORE PROPAGATING IT...
      // Don't bother propagating to z_s = 0 surface, just reject it now.
      //////////////////////////////////////////////////////////////////
      if ( GammaPhoton.m_p4.r0 < m_egy_gamma_min ) {
	//if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
	GammaPhoton.Zero();
	continue;
      }

      //////////////////////////////////////////////////////////////
      // ASK: Is energy of photon great enough to check and see if
      // it will cause a second generation of pair production?
      //
      // If so, then call PropagateSecondaryPhoton.
      // If not, call new function PropagatePhotonToObserver()
      //////////////////////////////////////////////////////////////
      double egy = Double(GammaPhoton.m_p4.r0);
      double redshift = Double(GammaPhoton.m_z_s);
      double tau = GetOpticalDepthVal(egy,redshift);
      if(tau < m_tauCutoff) {
	PropagatePhotonToObserver(GammaPhoton);
	if (GammaPhoton.m_p4.r0 >= m_egy_gamma_min)
	  StoreSecPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
			 GammaPhoton.m_weight);
	GammaPhoton.Zero();
	continue;
      } 
      
      //////Otherwise, just propagate IC scattered photon to z_s = 0////////
      RelParticle* Electron = new RelParticle;
      RelParticle* Positron = new RelParticle;
      
      bool pair_prod = 
	PropagateSecondaryPhoton(GammaPhoton,Electron,Positron);
      
      if (!pair_prod) {
	delete Electron;
	delete Positron;
	if(GammaPhoton.m_z_s > 1.0e-15) { // Should never happen
	  cerr<<"\nWARNING: Caught z_s=0 photon with z_s> 1E-15!"<<endl;
	  cerr<<"(Secondary photon...)"<<endl<<endl;
	}
	if (GammaPhoton.m_p4.r0 >= m_egy_gamma_min)
	  StoreSecPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
			 GammaPhoton.m_weight);
	if (m_trk_delay_bool) SaveToTrackTimeDelayFile(GammaPhoton);
	else {
	  //if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
	}
      }
      else { // Add electrons to stack?
	cout<<"---Secondary Photon Pair Production Triggered---"<<endl;
	if(Positron->m_p4.r0 >m_egy_lepton_min) {
	  m_globalLeptonNum++;
	  Positron->m_tag = m_globalLeptonNum;
	  lepton_stack.push(Positron);
	  if (m_trk_leptons_bool) {
	    CreateLeptonTree(Positron);
	    SaveLepton(Positron);
	  }
	}
	else {
	  //if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Positron);
	  delete Positron;
	}
	if(Electron->m_p4.r0>m_egy_lepton_min) {
	  m_globalLeptonNum++;
	  Electron->m_tag = m_globalLeptonNum;
	  lepton_stack.push(Electron);
	  if (m_trk_leptons_bool) {
	    CreateLeptonTree(Electron);
	    SaveLepton(Electron);
	  }
	}
	else {
	  //if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Electron);
	  delete Electron;
	}
      }
      if(lepton_stack.size() > stack_max) 
	stack_max = lepton_stack.size();
      GammaPhoton.Zero();
      //} // tau greater than tau cutoff, propagate secondary photon

	//} // photon is greater than E_min
      
    } // End while loop over leptons
      //PhotonList.close();
    cout<<"\nMax number of leptons in stack = "<<stack_max<<endl<<endl;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  void IGCascadeSim::CreateDirectGamma(RelParticle& GammaPhoton)
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


  void IGCascadeSim::
  PropagateDirectPhoton(RelParticle& GammaPhoton, 
			stack<RelParticle*>& lepton_stack)
  /*!  
    At the end of this function, the Electron and Positron are
    created at the interaction redshift, IF a pair production event
    takes place.
  */
  {

    double z_min = GetMinZ(GammaPhoton.m_p4, GammaPhoton.m_r4, GammaPhoton.m_z);
    double tot_lambda_int = 0.0;
    VEC3D_T z_int = "0.0";
    VEC3D_T z_emit = GammaPhoton.m_z;

    RelParticle* Electron = new RelParticle;
    RelParticle* Positron = new RelParticle;
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
      
    }
    else { 
      // z_min is not reached...However, from propagating photon
      // to z_final, it is possible that z_s=0 surface could be
      // crossed, since z_min does not correspond exactly to the
      // z_s=0 surface for computation time issues.
      VEC3D_T delta_z = (z_emit - z_int);
      bool pair_prod = 
	m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				    GammaPhoton.m_z,GammaPhoton.m_z_s, 
				    delta_z);
      if( pair_prod ) {
	// Pair production has legitimately occured and proceed as
	// usual defining EBL photon interacted with, Leptons, etc.
	RelParticle EBLPhoton;
	EBLPhoton.m_p4.r0 = 
	  m_pspace->GetEBLPhotonEgy(m_ebl,tot_lambda_int,
				    Double(z_int),Double(z_emit),
				    Double(GammaPhoton.m_p4.r0));
	m_pspace->UpdateEBLPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				  GammaPhoton.m_z,GammaPhoton.m_z_s,
				  EBLPhoton.m_p4,EBLPhoton.m_r4,EBLPhoton.m_z);
      
	////////////////////////////////////////////////////
	// CREATE AND DEFINE LEPTONS PRODUCED BY PAIRPROD //
	////////////////////////////////////////////////////
	DefineLeptons(GammaPhoton,EBLPhoton,Electron,Positron);
	GammaPhoton.Zero();

	// Decide whether or not to put the leptons onto the stack...
	if(Positron->m_p4.r0 > m_egy_lepton_min) {
	  m_globalLeptonNum++;
	  Positron->m_tag = m_globalLeptonNum;
	  lepton_stack.push(Positron);
	  if (m_trk_leptons_bool) {
	    CreateLeptonTree(Positron);
	    SaveLepton(Positron);
	  }
	}
	else {
	  //if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Positron);
	  delete Positron;
	}
	if(Electron->m_p4.r0 > m_egy_lepton_min) {
	  m_globalLeptonNum++;
	  Electron->m_tag = m_globalLeptonNum;
	  lepton_stack.push(Electron);
	  if (m_trk_leptons_bool) {
	    CreateLeptonTree(Electron);
	    SaveLepton(Electron);
	  }
	}
	else {
	  //if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Electron);
	  delete Electron;
	}
      }

    }
    

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGCascadeSim::
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
    std::cout<<"Pair Production Occurs..."<<std::endl;
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
    std::cout<<endl<<"new RelParticle: Electron..."<<endl;
    
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
    std::cout<<"new RelParticle: Positron..."<<endl<<endl;
    /////////////////////////////////////////

    //cout<<"gamma: "<<GammaPhoton.m_p4<<endl;
    //cout<<"EBLPhoton: "<<EBLPhoton.m_p4<<endl;
    
    bool RelKinematics = 
      m_pspace->RelativisticKinematics(GammaPhoton.m_p4,EBLPhoton.m_p4,
				       Electron->m_p4, Positron->m_p4);
    
    if (!RelKinematics)
      std::cerr<<std::endl<<"ERROR: RelativisticKinematics failed...\n\n";
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void IGCascadeSim::CreateLeptonTree(RelParticle* Lepton)
  {
    ostringstream oss_q;
    oss_q << Lepton->m_q;
    ostringstream oss_tag;
    oss_tag << Lepton->m_tag;
    cout<<"  tag: "<<Lepton->m_tag;
    string treeName = "electron_"+oss_tag.str()+"_q"+oss_q.str();
    cout<<"...Creating electron tree: "<<treeName<<endl;
    TTree* elecTree = new TTree(treeName.c_str(),"redshift,egy,px,py,pz,time,rx,ry,rz");
    //elecTree->SetDirectory(0);
    elecTree->Branch("redshift",&Lepton->m_tz,"redshift/D");
    elecTree->Branch("egy",&Lepton->m_tegy,"egy/D");
    elecTree->Branch("px",&Lepton->m_tpx,"px/D");
    elecTree->Branch("py",&Lepton->m_tpy,"py/D");
    elecTree->Branch("pz",&Lepton->m_tpz,"pz/D");
    elecTree->Branch("time",&Lepton->m_ttime,"time/D");
    elecTree->Branch("rx",&Lepton->m_trx,"rx/D");
    elecTree->Branch("ry",&Lepton->m_try,"ry/D");
    elecTree->Branch("rz",&Lepton->m_trz,"rz/D");

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void IGCascadeSim::WriteLeptonToFile(RelParticle* Lepton)
  {
    ostringstream oss_q;
    oss_q << Lepton->m_q;
    ostringstream oss_tag;
    oss_tag << Lepton->m_tag;
    string treeName = "electron_"+oss_tag.str()+"_q"+oss_q.str();
    TTree* elecTree = (TTree*)gDirectory->Get(treeName.c_str());
    TFile* rootfile = new TFile(m_save_lepton_file.c_str(),"UPDATE");
    elecTree->Write();
    rootfile->Close();
    delete rootfile;
    delete elecTree;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double IGCascadeSim::
  GetOpticalDepthVal(double egy, double z)
  {
    double TeV = 1.0e12;

    double egyMin = m_optDepthTable->GetRowVal(1);
    if(egy <= (egyMin*TeV)) return 0.0;
    int rowMax = m_optDepthTable->GetNRows();
    double egyMax = m_optDepthTable->GetRowVal(rowMax-1);
    if(egy >= (egyMax*TeV)) return 1.0;
    
    double zMin = m_optDepthTable->GetColVal(1);
    if(z <= zMin) return 1.0;
    int colMax = m_optDepthTable->GetNCols();
    double zMax = m_optDepthTable->GetColVal(colMax-1);
    if(z >= zMax) return 1.0;
    
    egy = egy/TeV;
    return m_optDepthTable->LinInterpolate(egy,z);
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  void IGCascadeSim::
  PropagatePhotonToObserver(RelParticle& GammaPhoton)
  {
    
    double z_min = GetMinZ(GammaPhoton.m_p4, GammaPhoton.m_r4, GammaPhoton.m_z);
    VEC3D_T z_emit = GammaPhoton.m_z;
    VEC3D_T delta_z = (z_emit - z_min);
    //bool pair_prod = 
    m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,GammaPhoton.m_r4,
				GammaPhoton.m_z, GammaPhoton.m_z_s, 
				delta_z);
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  bool IGCascadeSim::
  PropagateSecondaryPhoton(RelParticle& GammaPhoton, RelParticle* Electron, 
			   RelParticle* Positron)
  {

    //#ifdef DEBUG 
    //cout<<"Propagating SECONDARY photon through EBL...";
    //#endif
    if (m_single_gen_bool) {
      PropagatePhotonToObserver(GammaPhoton);
      return false;
    }
    double z_min = GetMinZ(GammaPhoton.m_p4, GammaPhoton.m_r4, GammaPhoton.m_z);
    VEC3D_T z_emit = GammaPhoton.m_z;

    //if (FORCE_SINGLE_GEN) {
    //VEC3D_T delta_z = (z_emit - z_min);
    //bool pair_prod = m_pspace->UpdateGammaPhoton(GammaPhoton.m_p4,
    //						   GammaPhoton.m_r4,
    //				  GammaPhoton.m_z,GammaPhoton.m_z_s, 
    //				  delta_z);
    //return false;
    //}
    
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

    //cout<<"Secondary Photon propagation complete...Parameters:\n";
    //cout<<"egy: "<<GammaPhoton.m_p4.r0<<" zs: "<<GammaPhoton.m_z_s<<endl;

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
    
    // CREATE AND DEFINE LEPTONS PRODUCED BY PAIRPROD
    DefineLeptons(GammaPhoton,EBLPhoton,Electron,Positron);     
    
    return true;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
  double IGCascadeSim::GetMinZ(const Vec4D& gam_ph_p4, const Vec4D& gam_ph_r4,
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


  void IGCascadeSim::PropagateLepton(stack<RelParticle*>& lepton_stack, 
				     RelParticle& GammaPhoton)
  {

    //#ifdef DEBUG
    //#endif

    RelParticle* Lepton = lepton_stack.top();
    lepton_stack.pop();
    //cout<<"Propagating Lepton, at z = "<<Lepton->m_z_s<<" egy: "<<
    //Lepton->m_p4.r0<<endl;
    //cout<<"Propagating Lepton, tag = "<<Lepton->m_tag<<" egy: "<<Lepton->m_p4.r0<<endl;
    
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
    
    if (m_trk_delay_bool) GammaPhoton.m_elec_time = Lepton->m_r4.r0;
    
    if(Lepton->m_p4.r0 > m_egy_lepton_min && 
       Lepton->m_z_s > m_BFieldGrid->m_DE) {
      lepton_stack.push(Lepton);
      if (m_trk_leptons_bool) SaveLepton(Lepton);
    } else {
      //if(SAVE_LOW_EGY) SaveToLowEnergyFile(*Lepton);
      cout<<"lepton z_s = "<<Lepton->m_z_s<<" egy: "<<Lepton->m_p4.r0
	  <<endl;
      if (m_trk_leptons_bool) WriteLeptonToFile(Lepton);
      cout<<"lepton "<<Lepton->m_tag<<" deleted from stack"<<endl<<endl;
      delete Lepton;
    }

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  void IGCascadeSim::
  SaveDirectPhoton(RelParticle& GammaPhoton, ofstream& PhotonList)
  {
        
    if(GammaPhoton.m_z_s > 1.0e-15) { // Should never happen
      cerr<<"\nWARNING: Caught z_s=0 photon with z_s > 1E-15!\n"<<endl;
      cerr<<"z_s: "<<GammaPhoton.m_z_s<<" egy: "<<GammaPhoton.m_p4.r0<<endl;
    }
    
    //if (GammaPhoton.m_p4.r0 >= m_egy_gamma_min) {
      //WritePhotonToFile(PhotonList, GammaPhoton.m_p4, 
      //		GammaPhoton.m_r4, GammaPhoton.m_weight);
    //}
    // else 
    //if(SAVE_LOW_EGY) SaveToLowEnergyFile(GammaPhoton);
    

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void IGCascadeSim::
  StoreSecPhoton(Vec4D& gam_ph_p4, Vec4D& gam_ph_r4,const double& weight)
  {

    //VEC3D_T EPS = 1.0E-25;
    Vec3D e_r = gam_ph_r4.r/gam_ph_r4.r.Norm();
    Vec3D e_p = gam_ph_p4.r/gam_ph_p4.r.Norm();
    m_theta   = 0.0;
    m_phi     = 0.0;
    VEC3D_T DE = "1.0e-25";
    // To account for the case of zero interactions:
    if(fabs(e_r.x) > m_DE || fabs(e_r.y) > m_DE ) { 
      m_theta = Double(atan2(sqrt(e_r.x*e_r.x + e_r.y*e_r.y),e_r.z));
      m_phi = Double(atan2(e_r.y,e_r.x));
    }
    m_thetap  = 0.0;
    double phi_p   = 0.0;
    if(fabs(e_p.x) > m_DE || fabs(e_p.y) > m_DE) {
      phi_p= Double(atan2(e_p.y,e_p.x));
      m_thetap = Double(atan2(sqrt(e_p.x*e_p.x + e_p.y*e_p.y),e_p.z));
    }
    m_xi = m_phi - phi_p;
    
    m_time = Double(gam_ph_r4.r0);
    m_egyPrim = Double(m_egy_cascade);
    m_egySec = Double(gam_ph_p4.r0);
    m_weight = weight;
    m_secPhotonTree->Fill();
    
    //stream<<gam_ph_p4.r0 <<" "<<theta<<" "<<phi<<" "<<
    //gam_ph_r4.r0<<" "<<theta_p<<" "<<xi<<" "<<weight <<std::endl;

    
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  void IGCascadeSim::SaveToLowEnergyFile(RelParticle& Particle)
  {
    
    ofstream low_egy_stream(m_low_egy_file.c_str(), std::ios::app);
    low_egy_stream<<Particle.m_p4.r0<<" "<<Particle.m_z<<" "<<Particle.m_q<<
      endl;
    low_egy_stream.close();
    // Nothing yet...    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  


  void IGCascadeSim::SaveLepton(RelParticle* Lepton)
  {
    
    // Saves the line:
    //   number_tag q z E px py pz time x y z
    //
    //ofstream ofile(m_save_lepton_file.c_str(), std::ios::app);
    //ofile<<Lepton->m_tag<<" "<<Lepton->m_q<<" "<<Lepton->m_z<<" "<<Lepton->m_p4.r0<<" "<<Lepton->m_p4.r<<" "<<Lepton->m_r4.r0<<" "<<Lepton->m_r4.r<<endl;
    //ofile.close();
    
    ostringstream oss_q;
    oss_q << Lepton->m_q;
    ostringstream oss_tag;
    oss_tag << Lepton->m_tag;
    string treeName = "electron_"+oss_tag.str()+"_q"+oss_q.str();
    TTree* elecTree = (TTree*)gDirectory->Get(treeName.c_str());
    //cout<<"  z: "<<Lepton->m_z<<" time: "<<Lepton->m_r4.r0<<endl; 
    //cout<<"  Filling tree: "<<treeName<<endl;
    Lepton->m_tz = Double(Lepton->m_z);
    Lepton->m_tegy = Double(Lepton->m_p4.r0);
    Lepton->m_tpx = Double(Lepton->m_p4.r.x);
    Lepton->m_tpy = Double(Lepton->m_p4.r.y);
    Lepton->m_tpz = Double(Lepton->m_p4.r.z);
    Lepton->m_ttime = Double(Lepton->m_r4.r0);
    Lepton->m_trx = Double(Lepton->m_r4.r.x);
    Lepton->m_try = Double(Lepton->m_r4.r.y);
    Lepton->m_trz = Double(Lepton->m_r4.r.z);
    //cout<<"--Saving Lepton..."<<Lepton->m_tegy<<" "<<Lepton->m_tz;
    //cout<<" "<<Lepton->m_tpx<<" "<<Lepton->m_tpy<<" "<<Lepton->m_tpz;
    //cout<<" "<<Lepton->m_ttime;
    //cout<<" "<<Lepton->m_trx<<" "<<Lepton->m_try<<" "<<Lepton->m_trz;
    
    elecTree->Fill();
    //cout<<"...Filled."<<endl;
    
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void IGCascadeSim::SaveToTrackTimeDelayFile(RelParticle& Photon)
  {
    ofstream ofile(m_track_time_delay_file.c_str(),std::ios::app);
    
    Vec3D e_r = Photon.m_r4.r/Photon.m_r4.r.Norm();
    Vec3D e_p = Photon.m_p4.r/Photon.m_p4.r.Norm();
    VEC3D_T theta   = "0.0";
    VEC3D_T phi     = "0.0";
    
    // To account for case of -> 0 angle
    if(fabs(e_r.x) > m_DE || fabs(e_r.y) > m_DE ) {
      theta = atan2(sqrt(e_r.x*e_r.x + e_r.y*e_r.y),e_r.z);
      phi = atan2(e_r.y,e_r.x);
    }
    
    VEC3D_T theta_p  = "0.0";
    VEC3D_T phi_p   = "0.0";
    if(fabs(e_p.x) > m_DE || fabs(e_p.y) > m_DE) {
      phi_p=atan2(e_p.y,e_p.x);
      theta_p = atan2(sqrt(e_p.x*e_p.x + e_p.y*e_p.y),e_p.z);
    }
    VEC3D_T xi = phi - phi_p;
    
    ofile<<std::fixed<<Photon.m_p4.r0 <<" "<<theta<<" "<<phi<<" "<< theta_p
	 <<" "<<xi<<" "<< Photon.m_r4.r0<<" "<< Photon.m_elec_time <<std::endl;
    
    ofile.close();
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
}

