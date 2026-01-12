//-*-mode:c++; mode:font-lock;-*-

#include<iostream>
#include<fstream>
#include<iterator>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<cctype>
#include<errno.h>
#include<cstring>

#include<EBLModelsFile.hpp>

union Magic
{
  unsigned i;
  char c[4];
};

const unsigned EBLModelsFile::scVersion = 1;

EBLModelsFile::~EBLModelsFile()
{
#ifdef LOUDLOUDLOUD
  std::cerr << "EBLModelsFile::~EBLModelsFile()" << std::endl;
#endif

  if(fWritable)
    {
      if(fseek(fFP, fSeekPos, SEEK_SET) < 0)
	{
	  perror("fseek [1]");
	  exit(EXIT_FAILURE);
	}

      if((fwrite(&fEBLNumModels, sizeof(fEBLNumModels), 1, fFP) != 1) ||
	 (fwrite(&fAGNNumModels, sizeof(fAGNNumModels), 1, fFP) != 1) ||
	 (fwrite(&fAcceptableModelsNum,sizeof(fAcceptableModelsNum),1,fFP)!=1))
	{
	  perror("fwrite [0]");
	  exit(EXIT_FAILURE);
	}
    }

  for(std::vector<AGNModel*>::iterator i = fAGNModels.begin();
      i != fAGNModels.end(); i++)delete *i;

  for(std::vector<ModelMatch*>::iterator i = fAcceptableModels.begin();
      i != fAcceptableModels.end(); i++)delete *i;

  fclose(fFP);

  delete[] fEBLParameters;
}

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// W   W  RRRR   III  TTTTT  EEEEE
// W   W  R   R   I     T    E
// W W W  RRRR    I     T    EEE
// W W W  R   R   I     T    E
//  W W   R   R  III    T    EEEEE

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

EBLModelsFile::EBLModelsFile(const std::string& filename,
			     unsigned num_ebl_param, const double* ebl_lambda,
			     unsigned num_agn_param,
			     const char** agn_param_names,
			     unsigned num_match_param,
			     const char** match_param_names,
			     const std::string& vhe_dataset_name,
			     unsigned vhe_num_data,
			     const double* vhe_e,
			     const double* vhe_i, const double* vhe_di):
  fFileName(filename), fFP(), fWritable(true), fSeekPos(0),
  fEBLNumParameters(num_ebl_param), fEBLLambda(), fEBLParameters(0),
  fEBLTau01(), fEBLTau1(), fEBLTau10(), fEBLIStars(), fEBLIDust(),
  fEBLNumModels(0),
  fAGNNumParameters(num_agn_param), fAGNParameterNames(),
  fAGNModelsLookup(),  fAGNModels(), fAGNModelsUnwritten(),
  fAGNNumModels(0),
  fAcceptableNumParameters(num_match_param), fAcceptableParameterNames(),
  fAcceptableModels(), fNextAcceptableModel(), fAcceptableModelsNum(0),
  fVHEDatasetName(vhe_dataset_name), fVHENumData(vhe_num_data),
  fVHEDataE(), fVHEDataI(), fVHEDataIError()
{
  // **************************************************************************
  // Set the class member variables
  // **************************************************************************

  for(unsigned i=0;i<num_ebl_param;i++)
    fEBLLambda.push_back(ebl_lambda[i]);

  for(unsigned i=0;i<num_agn_param;i++)
    fAGNParameterNames.push_back(std::string(agn_param_names[i]));

  for(unsigned i=0;i<num_match_param;i++)
    fAcceptableParameterNames.push_back(std::string(match_param_names[i]));

  for(unsigned i=0;i<vhe_num_data;i++)
    {
      fVHEDataE.push_back(vhe_e[i]);
      fVHEDataI.push_back(vhe_i[i]);
      fVHEDataIError.push_back(vhe_di[i]);
    }

  fEBLParameters = new double[num_ebl_param];

  // **************************************************************************
  // Open file and write the data
  // **************************************************************************

  fFP = fopen(filename.c_str(), "w");
  if(fFP == NULL)
    {
      std::cerr << "Could not open file " << filename << std::endl
		<< "fopen: " << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 1) Write the MAGIC number
  // --------------------------------------------------------------------------
  Magic magic;
  magic.c[0] = 'E';
  magic.c[1] = 'B';
  magic.c[2] = 'L';
  magic.c[3] = '\0';
  if(!b_write(magic.i, fFP))
    {
      perror("write [1]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 2) Write the VERSION number
  // --------------------------------------------------------------------------
  if(!b_write(scVersion, fFP))
    {
      perror("write [2]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 3) Write the VHE dataset name
  // --------------------------------------------------------------------------
  if(!b_write(vhe_dataset_name, fFP))
    {
      perror("write [3]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 4) Write the VHE dataset
  // --------------------------------------------------------------------------
  if(!b_write(fVHENumData, fFP) ||
     !b_write(fVHEDataE, fFP) ||
     !b_write(fVHEDataI, fFP) ||
     !b_write(fVHEDataIError, fFP))
    {
      perror("write [4]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 5) Write the number of EBL, AGN and MATCH parameters
  // --------------------------------------------------------------------------
  if(!b_write(fEBLNumParameters, fFP) ||
     !b_write(fAGNNumParameters, fFP) ||
     !b_write(fAcceptableNumParameters, fFP))
    {
      perror("write [5]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 6) Write the EBL wavelengths and AGN and MATCH parameter names
  // --------------------------------------------------------------------------
  if(!b_write(fEBLLambda, fFP))
    {
      perror("write [6]");
      exit(EXIT_FAILURE);
    }

  if(!b_write(fAGNParameterNames, fFP))
    {
      perror("write [7]");
      exit(EXIT_FAILURE);
    }

  if(!b_write(fAcceptableParameterNames, fFP))
    {
      perror("write [8]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 7) Write the number of EBL and AGN models
  // --------------------------------------------------------------------------

  // We write zeros here now and seek back at the end to write real values
  fSeekPos = ftell(fFP);
  if(!b_write(fEBLNumModels, fFP) ||
     !b_write(fAGNNumModels, fFP) ||
     !b_write(fAcceptableModelsNum, fFP))
    {
      perror("write [9]");
      exit(EXIT_FAILURE);
    }
}

void EBLModelsFile::setEBLModel(const double* param,
				double tau01, double tau1, double tau10,
				double i_stars, double i_dust)
{
  assert(fWritable);

  std::copy(param, param+fEBLNumParameters, fEBLParameters);
  fEBLTau01  = tau01;
  fEBLTau1   = tau1;
  fEBLTau10  = tau10;
  fEBLIStars = i_stars;
  fEBLIDust  = i_dust;
  fEBLNumModels++;

  for(std::vector<ModelMatch*>::iterator i = fAcceptableModels.begin();
      i != fAcceptableModels.end(); i++)delete *i;
  fAcceptableModels.clear();
}

void EBLModelsFile::writeEBLModel()
{
#if 0
  std::cout
    << "EBL models written:       " << fEBLNumModels << std::endl
    << "AGN Models registered:    " << fAGNNumModels << std::endl
    << "AGN parameters size:      " << fAGNModels.size() << std::endl
    << "AGN parameters unwritten: " << fAGNModelsUnwritten.size() << std::endl
    << "AGN lookup size:          " << fAGNModelsLookup.size() << std::endl;

  unsigned count = 0;
  for(std::set<AGNModel*,AGNModel>::iterator i =
	fAGNModelsLookup.begin(); i != fAGNModelsLookup.end(); i++)
    count++;

  std::cout
    << "AGN lookup count:         " << count << std::endl;
#endif

  // --------------------------------------------------------------------------
  // Check consistancy
  // --------------------------------------------------------------------------
  assert(fWritable);

  // --------------------------------------------------------------------------
  // 1) Write the AGN models accumulated
  // --------------------------------------------------------------------------

  if(!b_write(fAGNModelsUnwritten, fFP))
    {
      perror("write [10]");
      exit(EXIT_FAILURE);
    }

  fAGNModelsUnwritten.clear();

  // --------------------------------------------------------------------------
  // 2) Write the EBL parameters
  // --------------------------------------------------------------------------

  if(!b_write(fEBLParameters, fEBLNumParameters, fFP) ||
     !b_write(fEBLTau01, fFP) ||
     !b_write(fEBLTau1, fFP) ||
     !b_write(fEBLTau10, fFP) ||
     !b_write(fEBLIStars, fFP) ||
     !b_write(fEBLIDust, fFP))
    {
      perror("write [11]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 3) Write the acceptable models
  // --------------------------------------------------------------------------
  if(!b_write(fAcceptableModels, fFP))
    {
      perror("write [12]");
      exit(EXIT_FAILURE);
    }
}

void EBLModelsFile::setAcceptableAGNModel(const double* param,
					  const float* match_param)
{
  assert(fWritable);

  // Find the model using the lookup table (std::set -- a binary search tree)
  AGNModel key(0,fAGNNumParameters,param);
  AGNModel* model = 0;
  std::set<AGNModel*,AGNModel>::iterator found = fAGNModelsLookup.find(&key);

  if(found == fAGNModelsLookup.end())
    {
      // Model not found -- so add it now
      model = new AGNModel(fAGNNumModels++, fAGNNumParameters, param);
      fAGNModels.push_back(model);
      fAGNModelsLookup.insert(model);

#if 0
      if(fAGNModelsLookup.size() != fAGNModels.size())
	{
	  std::cout << "Found mismatch for parameters.." << std::endl;
	  for(unsigned i=0;i<fAGNNumParameters;i++)
	    std::cout << param[i] << ' ' << model->fParam[i] << std::endl;

	  std::cout << "Dumping fAGNModels: " << std::endl;
	  std::copy(fAGNModels.begin(), fAGNModels.end(),
		    std::ostream_iterator<const AGNModel*>(std::cout, "\n"));

	  std::cout << "Dumping fAGNModelsLookup: " << std::endl;
	  std::copy(fAGNModelsLookup.begin(), fAGNModelsLookup.end(),
		    std::ostream_iterator<const AGNModel*>(std::cout, "\n"));

	  assert(0);
	}
#endif

      fAGNModelsUnwritten.push_back(model);
    }
  else
    {
      // Model found
      model = *found;
    }

#if 0
  std::cerr << fEBLNumModels << '\t' << fAGNNextIndex << '\t';
  for(unsigned i=0;i<fAGNNumParameters;i++)
    std::cerr << param[i] << '\t' << fAGNModels[fAGNNextIndex][i] << '\t';
  std::cerr << std::endl;
#endif

  fAcceptableModels.push_back(new ModelMatch(model->fModelIndex,
					     fAcceptableNumParameters,
					     match_param));
  fAcceptableModelsNum++;
}


// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// RRRR   EEEEE   AAA   DDDD
// R   R  E      A   A  D   D
// RRRR   EEE    AAAAA  D   D
// R   R  E      A   A  D   D
// R   R  EEEEE  A   A  DDDD

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

EBLModelsFile::EBLModelsFile(const std::string& filename):
  fFileName(filename), fFP(), fWritable(false), fSeekPos(0),
  fEBLNumParameters(0), fEBLLambda(), fEBLParameters(0),
  fEBLTau01(), fEBLTau1(), fEBLTau10(), fEBLIStars(), fEBLIDust(),
  fEBLNumModels(0),
  fAGNNumParameters(), fAGNParameterNames(),
  fAGNModelsLookup(),  fAGNModels(), fAGNModelsUnwritten(),
  fAGNNumModels(0),
  fAcceptableModels(), fNextAcceptableModel(),
  fVHEDatasetName(), fVHENumData(),
  fVHEDataE(), fVHEDataI(), fVHEDataIError()
{
  // **************************************************************************
  // Open file and read the header data
  // **************************************************************************

  fFP = fopen(filename.c_str(), "r");
  if(fFP == NULL)
    {
      std::cerr << "Could not open file " << filename << std::endl
		<< "fopen: " << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 1) Read and test the MAGIC number
  // --------------------------------------------------------------------------
  Magic magic;
  if(!b_read(magic.i, fFP))
    {
      perror("read [1]");
      exit(EXIT_FAILURE);
    }

  assert((magic.c[0] == 'E') && (magic.c[1] == 'B') && (magic.c[2] == 'L') &&
	 (magic.c[3] == '\0'));

  // --------------------------------------------------------------------------
  // 2) Read and test the VERSION number
  // --------------------------------------------------------------------------
  unsigned version;
  if(!b_read(version, fFP))
    {
      perror("read [2]");
      exit(EXIT_FAILURE);
    }

  assert(version == scVersion);

  // --------------------------------------------------------------------------
  // 3) Read the VHE dataset name
  // --------------------------------------------------------------------------
  if(!b_read(fVHEDatasetName, fFP))
    {
      perror("read [3]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 4) Read the VHE dataset
  // --------------------------------------------------------------------------
  if(!b_read(fVHENumData, fFP) ||
     !b_read(fVHEDataE, fFP) ||
     !b_read(fVHEDataI, fFP) ||
     !b_read(fVHEDataIError, fFP))
    {
      perror("read [4]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 5) Read the number of EBL, AGN and MATCH parameters
  // --------------------------------------------------------------------------
  if(!b_read(fEBLNumParameters, fFP) ||
     !b_read(fAGNNumParameters, fFP) ||
     !b_read(fAcceptableNumParameters, fFP))
    {
      perror("read [5]");
      exit(EXIT_FAILURE);
    }

  fEBLParameters = new double[fEBLNumParameters];

  // --------------------------------------------------------------------------
  // 6) Read the EBL wavelengths and AGN and MATCH parameter names
  // --------------------------------------------------------------------------
  if(!b_read(fEBLLambda, fFP))
    {
      perror("read [6]");
      exit(EXIT_FAILURE);
    }

  if(!b_read(fAGNParameterNames, fFP))
    {
      perror("read [7]");
      exit(EXIT_FAILURE);
    }

  if(!b_read(fAcceptableParameterNames, fFP))
    {
      perror("read [8]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 7) Read the number of EBL and AGN models
  // --------------------------------------------------------------------------

  if(!b_read(fEBLNumModels, fFP) ||
     !b_read(fAGNNumModels, fFP) ||
     !b_read(fAcceptableModelsNum, fFP))
    {
      perror("read [9]");
      exit(EXIT_FAILURE);
    }

#if 0
  std::cerr << fEBLNumModels << std::endl
	    << fAGNNumModels << std::endl
	    << fAcceptableModelsNum << std::endl;
#endif

  // --------------------------------------------------------------------------
  // 8) Read through the file and extract all AGN models
  // --------------------------------------------------------------------------

  // Remember where to return later
  fSeekPos = ftell(fFP);

  unsigned num_ebl_models = 0;
  unsigned num_acceptable_models = 0;

  while(1)
    {
      if(!b_read(fAGNModelsUnwritten, fAGNNumParameters, fFP))
	{
	  if(feof(fFP) && (fAGNModelsUnwritten.size() == 0))break;
	  perror("fread [10]");
	  exit(EXIT_FAILURE);
	}

      for(std::vector<AGNModel*>::iterator i = fAGNModelsUnwritten.begin();
	  i != fAGNModelsUnwritten.end(); i++)
	{
#if 0
	  std::cout << *i << std::endl;
#endif
	  AGNModel* model =
	    new AGNModel((*i)->fModelIndex, (*i)->fNumParam, (*i)->fParam);
	  fAGNModels.push_back(model);
	  fAGNModelsLookup.insert(model);
	}

      // SKIP the EBL model information
      if(fseek(fFP,
	       sizeof(double)*fEBLNumParameters+
	       sizeof(fEBLTau01)+sizeof(fEBLTau1)+sizeof(fEBLTau10)+
	       sizeof(fEBLIStars)+sizeof(fEBLIDust),
	       SEEK_CUR) < 0)
	{
	  perror("fseek [2]");
	  exit(EXIT_FAILURE);
	}

      // READ the number of matches and SKIP them
      unsigned num_acceptable;
      if(!b_read(num_acceptable, fFP))
	{
	  perror("read [11]");
	  exit(EXIT_FAILURE);
	}

      num_acceptable_models += num_acceptable;

      if(fseek(fFP,
	       num_acceptable*(sizeof(unsigned) +
			       fAcceptableNumParameters*sizeof(float)),
	       SEEK_CUR) < 0)
	{
	  perror("fseek [3]");
	  exit(EXIT_FAILURE);
	}

      num_ebl_models++;
    }

  // --------------------------------------------------------------------------
  // Sanity Checks
  // --------------------------------------------------------------------------

  assert(num_ebl_models == fEBLNumModels);
  assert(fAGNModels.size() == fAGNNumModels);
  assert(num_acceptable_models == fAcceptableModelsNum);

  // --------------------------------------------------------------------------
  // Rewind the file
  // --------------------------------------------------------------------------

  if(fseek(fFP, fSeekPos, SEEK_SET) < 0)
    {
      perror("fseek [4]");
      exit(EXIT_FAILURE);
    }

#if 0
  std::copy(fAcceptableParameterNames.begin(), fAcceptableParameterNames.end(),
	    std::ostream_iterator<std::string>(std::cout, "\n"));
#endif
}

bool EBLModelsFile::readNextEBLModel()
{
  // --------------------------------------------------------------------------
  // 1) SKIP over AGN models records
  // --------------------------------------------------------------------------

  unsigned num_agn_to_skip;

  bool temp = b_read(num_agn_to_skip, fFP);
  if(!temp)
    {
      if(feof(fFP))return false;
      else
	{
	  perror("fread [12]");
	  exit(EXIT_FAILURE);
	}
    }

  if(fseek(fFP,
	   num_agn_to_skip*(sizeof(unsigned)+sizeof(double)*fAGNNumParameters),
	   SEEK_CUR) < 0)
    {
      perror("fseek [5]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 2) Read the EBL parameters
  // --------------------------------------------------------------------------

  if(!b_read(fEBLParameters, fEBLNumParameters, fFP) ||
     !b_read(fEBLTau01, fFP) ||
     !b_read(fEBLTau1, fFP) ||
     !b_read(fEBLTau10, fFP) ||
     !b_read(fEBLIStars, fFP) ||
     !b_read(fEBLIDust, fFP))
    {
      perror("read [13]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // 3) Read the acceptable models
  // --------------------------------------------------------------------------

  if(!b_read(fAcceptableModels, fAcceptableNumParameters, fFP))
    {
      perror("read [14]");
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // Set iterator
  // --------------------------------------------------------------------------

  fNextAcceptableModel = fAcceptableModels.begin();
  return true;
}

bool EBLModelsFile::
nextAcceptableAGNModel(unsigned& param_index, const double* &param,
		       const float* &match_param)
{
  if(fNextAcceptableModel == fAcceptableModels.end())return false;
  param_index = (*fNextAcceptableModel)->fModelIndex;
  assert(param_index < fAGNNumModels);
  param = fAGNModels[param_index]->fParam;
  match_param = (*fNextAcceptableModel)->fParam;
  fNextAcceptableModel++;
  return true;
}

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// VHE DATA READER

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

VHEData::VHEData(const std::string& filename, double deltae_over_e):
  fEMultiplier(1.0+deltae_over_e), fName(filename), fParameters(),
  fEnergy(), fE2dNdE(), fDeltaE2dNdE(), fCommentToParameters(true)
{
  std::ifstream stream(filename.c_str());
  if(!stream)
    {
      perror("ifstream");
      exit(EXIT_FAILURE);
    }

  std::string line;
  while(getline(stream,line))
    {
      double e;
      double i;
      double di;
      std::istringstream linestream(line);
      linestream >> e >> i >> di;
      if(linestream)
	{
	  fEnergy.push_back(e*fEMultiplier);
	  fE2dNdE.push_back(i);
	  fDeltaE2dNdE.push_back(di);
	}
      else
	{
	  std::cerr << "VHEData::VHEData: parsing error in file " << filename
		    << std::endl
		    << line << std::endl;
	}
    }
}

bool VHEData::getline(std::istream& stream, std::string& line)
{
  line.clear();

  do
    {
      static const char* SPC = " \t\r\n";

      char c;
      stream.get(c); while((stream)&&(c!='\n')){ line+=c; stream.get(c); }
      if(line.empty())continue;

      // SKIP LEADING AND TRAILING SPACE
      std::string::size_type first_char = line.find_first_not_of(SPC);
      std::string::size_type last_char = line.find_last_not_of(SPC);
      line = line.substr(first_char, last_char+1);
      if(line.empty())continue;

      // SEPARATE COMMENT
      std::string comment;
      std::string::size_type hash = line.find('#');
      if(hash != std::string::npos)
	{
	  comment = line.substr(hash+1);
	  line = line.substr(0,hash);

	  if(!line.empty())
	    {
	      last_char = line.find_last_not_of(SPC);
	      line = line.substr(0, last_char+1);
	    }

	  if(!comment.empty())
	    {
	      first_char = comment.find_first_not_of(SPC);
	      comment = comment.substr(first_char);
	    }
	}

      // EXTRACT "#KEY:VALUE" DATA FROM TOP OF FILE
      if((line.empty()) && (!comment.empty()) && (fCommentToParameters))
	{
	  hash = comment.find('#');
	  if(hash != std::string::npos)
	    {
	      comment = comment.substr(0,hash);
	      if(!comment.empty())
		{
		  last_char = comment.find_last_not_of(SPC);
		  comment = comment.substr(0, last_char+1);
		}
	    }

	  std::string::size_type colon = comment.find(':');
	  if(colon != std::string::npos)
	    {
	      std::string key = comment.substr(0,colon);
	      std::string val = comment.substr(colon+1);

	      if(!key.empty())
		{
		  last_char = key.find_last_not_of(SPC);
		  key = key.substr(0,++last_char);
		  for(std::string::iterator i=key.begin(); i!=key.end(); i++)
		    *i = char(tolower(*i));
		}

	      if(!val.empty())
		{
		  first_char = val.find_first_not_of(SPC);
		  val = val.substr(first_char);
		}

	      if(!key.empty())fParameters[key]=val;
	    }
	  else fCommentToParameters = false;
	}
      else fCommentToParameters = false;
    }while((stream)&&(line.empty()));

  return !line.empty();
}

bool VHEData::redshift(double& z) const
{
  std::string z_string;
  if(!parameter("redshift",z_string))return false;
  std::istringstream z_stream(z_string);
  if(z_stream >> z)return true;
  else return false;
}

std::string VHEData::name() const
{
  std::ostringstream stream;
  stream << fName << " [x" << fEMultiplier << ']';
  return stream.str();
}

bool VHEData::parameter(const std::string& key, std::string& val) const
{

  std::map<std::string, std::string>::const_iterator i = fParameters.find(key);
  if(i == fParameters.end())return false;
  val = (*i).second;
  return true;
}
