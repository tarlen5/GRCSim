//-*-mode:c++; mode:font-lock;-*-

#ifndef EBLMODELSFILE_HPP
#define EBLMODELSFILE_HPP

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

template <class T> struct Model {
  unsigned fModelIndex;
  unsigned fNumParam;
  T *fParam;
  Model() : fModelIndex(0), fNumParam(0), fParam(0) {
#ifdef LOUDLOUDLOUD
    std::cerr << "CREATE: " << (void *)this << ' ' << *this << ' ' << fParam
              << std::endl;
#endif
  }

  Model(unsigned idx, unsigned nparam, const T *param = 0)
      : fModelIndex(idx), fNumParam(nparam), fParam(0) {
    if (nparam) {
      fParam = new T[nparam];
      if (param) memcpy(fParam, param, nparam * sizeof(T));
    }
#ifdef LOUDLOUDLOUD
    std::cerr << "CREATE: " << (void *)this << ' ' << *this << ' ' << fParam
              << std::endl;
#endif
  }

  Model(const Model &o)
      : fModelIndex(o.fModelIndex), fNumParam(o.fNumParam), fParam(0) {
    if ((o.fParam) && (fNumParam)) {
      fParam = new T[fNumParam];
      memcpy(fParam, o.fParam, fNumParam * sizeof(T));
    }
#ifdef LOUDLOUDLOUD
    std::cerr << "CREATE: " << (void *)this << ' ' << *this << ' ' << fParam
              << std::endl;
#endif
  }

  const Model &operator=(const Model &o) {
    fModelIndex = o.fModelIndex;
    fNumParam = o.fNumParam;
    fParam = 0;
    if (o.fParam) {
      fParam = new T[fNumParam];
      memcpy(fParam, o.fParam, fNumParam * sizeof(T));
    }
#ifdef LOUDLOUDLOUD
    std::cerr << "EQUALS: " << (void *)this << ' ' << *this << ' ' << fParam
              << std::endl;
#endif
    return *this;
  }

  ~Model() {
#ifdef LOUDLOUDLOUD
    std::cerr << "DELETE: " << (void *)this << ' ' << *this << ' ' << fParam
              << std::endl;
#endif
    if (fParam) delete[] fParam;
  }

  bool operator<(const Model &o) const {
    assert(fNumParam == o.fNumParam);
    for (unsigned i = 0; i < fNumParam; i++)
      if (fParam[i] < (o.fParam[i] - 0.001))
        return true;
      else if (fParam[i] > (o.fParam[i] - 0.001))
        return false;
    return false;
  }

  bool operator()(const Model *a, const Model *b) const { return *a < *b; }

private:
};

template <class T>
inline std::ostream &operator<<(std::ostream &stream, const Model<T> &x) {
  stream << x.fModelIndex << ' ' << x.fNumParam;
  for (unsigned i = 0; i < x.fNumParam; i++)
    stream << ' ' << x.fParam[i];
  return stream;
}

template <class T>
inline std::ostream &operator<<(std::ostream &stream, const Model<T> *x) {
  return (stream << *x);
}

typedef Model<double> AGNModel;
typedef Model<float> ModelMatch;

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// EBLMODELSFILE

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

class EBLModelsFile {
public:
  virtual ~EBLModelsFile();

  // **************************************************************************
  // OUTPUT member functions
  // **************************************************************************

  EBLModelsFile(const std::string &filename, unsigned num_ebl_param,
                const double *ebl_lambda, unsigned num_agn_param,
                const char **agn_param_names, unsigned num_match_param,
                const char **match_param_names,
                const std::string &vhe_dataset_name, unsigned vhe_num_data,
                const double *vhe_e, const double *vhe_i, const double *vhe_di);

  void setEBLModel(const double *param, double tau01, double tau1, double tau10,
                   double i_stars, double i_dust);

  void writeEBLModel();

  void setAcceptableAGNModel(const double *param, const float *match_param);

  // **************************************************************************
  // INPUT member functions
  // **************************************************************************

  EBLModelsFile(const std::string &filename);
  bool readNextEBLModel();
  bool nextAcceptableAGNModel(unsigned &param_index, const double *&param,
                              const float *&match_param);

  // **************************************************************************
  // ACCESSORS
  // **************************************************************************

  unsigned eblNumParameters() const { return fEBLNumParameters; }
  const std::vector<double> &eblLambda() const { return fEBLLambda; }
  const double *eblParameters() const { return fEBLParameters; }
  double eblTau01() const { return fEBLTau01; }
  double eblTau1() const { return fEBLTau1; }
  double eblTau10() const { return fEBLTau10; }
  double eblIStars() const { return fEBLIStars; }
  double eblIDust() const { return fEBLIDust; }
  unsigned eblNumModels() const { return fEBLNumModels; }

  unsigned agnNumParameters() const { return fAGNNumParameters; }
  const std::vector<std::string> &agnParameterNames() const {
    return fAGNParameterNames;
  }
  unsigned agnNumModels() const { return fAGNModels.size(); }
  const double *agnModel(unsigned idx) const { return fAGNModels[idx]->fParam; }

  const std::string &vheDatasetName() const { return fVHEDatasetName; }
  unsigned vheNumData() const { return fVHENumData; }
  const std::vector<double> &vheDataE() const { return fVHEDataE; }
  const std::vector<double> &vheDataI() const { return fVHEDataI; }
  const std::vector<double> &vheDataIError() const { return fVHEDataIError; }

  unsigned acceptableNumParameters() const { return fAcceptableNumParameters; }
  const std::vector<std::string> &acceptableParameterNames() const {
    return fAcceptableParameterNames;
  }
  unsigned acceptableNumModels() const { return fAcceptableModelsNum; }

private:
  std::string fFileName;
  FILE *fFP;
  bool fWritable;
  long fSeekPos;

  unsigned fEBLNumParameters;
  std::vector<double> fEBLLambda;
  double *fEBLParameters;
  double fEBLTau01;
  double fEBLTau1;
  double fEBLTau10;
  double fEBLIStars;
  double fEBLIDust;
  unsigned fEBLNumModels;

  unsigned fAGNNumParameters;
  std::vector<std::string> fAGNParameterNames;
  std::set<AGNModel *, AGNModel> fAGNModelsLookup;
  std::vector<AGNModel *> fAGNModels;
  std::vector<AGNModel *> fAGNModelsUnwritten;
  unsigned fAGNNumModels;

  unsigned fAcceptableNumParameters;
  std::vector<std::string> fAcceptableParameterNames;
  std::vector<ModelMatch *> fAcceptableModels;
  std::vector<ModelMatch *>::const_iterator fNextAcceptableModel;
  unsigned fAcceptableModelsNum;

  std::string fVHEDatasetName;
  unsigned fVHENumData;
  std::vector<double> fVHEDataE;
  std::vector<double> fVHEDataI;
  std::vector<double> fVHEDataIError;

  static const unsigned scVersion;
};

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// BINARY INPUT-OUTPUT UTILITY FUNCTIONS

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// Should use a template here but VVV would not be happy!
inline bool b_write(const unsigned &o, FILE *fp);
inline bool b_write(const float &o, FILE *fp);
inline bool b_write(const double &o, FILE *fp);
inline bool b_write(const std::string &o, FILE *fp);
inline bool b_write(const char *o, FILE *fp);
inline bool b_write(const float *o, unsigned count, FILE *fp);
inline bool b_write(const double *o, unsigned count, FILE *fp);
inline bool b_write(const std::vector<double> &o, FILE *fp);
inline bool b_write(const std::vector<std::string> &o, FILE *fp);
inline bool b_write(const std::vector<AGNModel *> &o, FILE *fp);
inline bool b_write(const std::vector<ModelMatch *> &o, FILE *fp);

inline bool b_read(unsigned &o, FILE *fp);
inline bool b_read(float &o, FILE *fp);
inline bool b_read(double &o, FILE *fp);
inline bool b_read(std::string &o, FILE *fp);
// inline bool b_read(char* o, FILE* fp);
inline bool b_read(float *o, unsigned count, FILE *fp);
inline bool b_read(double *o, unsigned count, FILE *fp);
inline bool b_read(std::vector<double> &o, FILE *fp);
inline bool b_read(std::vector<std::string> &o, FILE *fp);
inline bool b_read(std::vector<AGNModel *> &o, unsigned num_param, FILE *fp);
inline bool b_read(std::vector<ModelMatch *> &o, unsigned num_param, FILE *fp);

inline bool b_write(const unsigned &o, FILE *fp) {
  return fwrite(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_write(const float &o, FILE *fp) {
  return fwrite(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_write(const double &o, FILE *fp) {
  return fwrite(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_write(const std::string &o, FILE *fp) {
  if (!b_write(unsigned(o.length()), fp)) return false;
  if (o.length() == 0) return true;
  return fwrite(o.data(), sizeof(*o.data()), o.length(), fp) == o.length();
}

inline bool b_write(const char *o, FILE *fp) {
  return b_write(std::string(o), fp);
}

inline bool b_write(const float *o, unsigned count, FILE *fp) {
  if (count == 0) return true;
  return fwrite(o, sizeof(*o), count, fp) == count;
}

inline bool b_write(const double *o, unsigned count, FILE *fp) {
  if (count == 0) return true;
  return fwrite(o, sizeof(*o), count, fp) == count;
}

inline bool b_write(const std::vector<double> &o, FILE *fp) {
  if (!b_write(unsigned(o.size()), fp)) return false;
  return b_write(&(*o.begin()), o.size(), fp);
}

inline bool b_write(const std::vector<std::string> &o, FILE *fp) {
  if (!b_write(unsigned(o.size()), fp)) return false;
  for (std::vector<std::string>::const_iterator i = o.begin(); i != o.end();
       i++)
    if (!b_write(*i, fp)) return false;
  return true;
}

inline bool b_write(const std::vector<AGNModel *> &o, FILE *fp) {
  if (!b_write(unsigned(o.size()), fp)) return false;
  for (std::vector<AGNModel *>::const_iterator i = o.begin(); i != o.end();
       i++) {
    if (!b_write((*i)->fModelIndex, fp)) return false;
    if (!b_write((*i)->fParam, (*i)->fNumParam, fp)) return false;
  }
  return true;
}

inline bool b_write(const std::vector<ModelMatch *> &o, FILE *fp) {
  if (!b_write(unsigned(o.size()), fp)) return false;
  if (o.size() == 0) return true;
  for (std::vector<ModelMatch *>::const_iterator i = o.begin(); i != o.end();
       i++) {
    if (!b_write((*i)->fModelIndex, fp)) return false;
    if (!b_write((*i)->fParam, (*i)->fNumParam, fp)) return false;
  }
  return true;
}

inline bool b_read(unsigned &o, FILE *fp) {
  return fread(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_read(float &o, FILE *fp) {
  return fread(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_read(double &o, FILE *fp) {
  return fread(&o, sizeof(o), 1, fp) == 1;
}

inline bool b_read(std::string &o, FILE *fp) {
  unsigned s;
  o.erase();
  if (!b_read(s, fp)) return false;
  if (s == 0) return true;
  char *buffer = new char[s];
  if (fread(buffer, sizeof(*buffer), s, fp) != s) {
    delete[] buffer;
    return false;
  }
  o.append(buffer, s);
  delete[] buffer;
  return true;
}

inline bool b_read(float *o, unsigned count, FILE *fp) {
  if (count == 0) return true;
  return fread(o, sizeof(*o), count, fp) == count;
}

inline bool b_read(double *o, unsigned count, FILE *fp) {
  if (count == 0) return true;
  return fread(o, sizeof(*o), count, fp) == count;
}

inline bool b_read(std::vector<double> &o, FILE *fp) {
  unsigned s;
  o.clear();
  if (!b_read(s, fp)) return false;
  if (s == 0) return true;
  double *buffer = new double[s];
  if (fread(buffer, sizeof(*buffer), s, fp) != s) {
    delete[] buffer;
    return false;
  }
  o.insert(o.begin(), buffer, buffer + s);
  delete[] buffer;
  return true;
}

inline bool b_read(std::vector<std::string> &o, FILE *fp) {
  unsigned s;
  o.clear();
  if (!b_read(s, fp)) return false;
  if (s == 0) return true;
  o.resize(s);
  for (std::vector<std::string>::iterator i = o.begin(); i != o.end(); i++) {
    if (!b_read(*i, fp)) return false;
  }
  return true;
}

inline bool b_read(std::vector<AGNModel *> &o, unsigned num_param, FILE *fp) {
  unsigned s;
  for (std::vector<AGNModel *>::iterator i = o.begin(); i != o.end(); i++)
    delete *i;
  o.clear();
  if (!b_read(s, fp)) return false;
  for (unsigned i = 0; i < s; i++) {
    unsigned index;
    if (!b_read(index, fp)) return false;
    AGNModel *model = new AGNModel(index, num_param, 0);
    if (!b_read(model->fParam, model->fNumParam, fp)) return false;
    o.push_back(model);
  }
  return true;
}

inline bool b_read(std::vector<ModelMatch *> &o, unsigned num_param, FILE *fp) {
  unsigned s;
  for (std::vector<ModelMatch *>::iterator i = o.begin(); i != o.end(); i++)
    delete *i;
  o.clear();
  if (!b_read(s, fp)) return false;
  for (unsigned i = 0; i < s; i++) {
    unsigned index;
    if (!b_read(index, fp)) return false;
    ModelMatch *model = new ModelMatch(index, num_param, 0);
    if (!b_read(model->fParam, model->fNumParam, fp)) return false;
    o.push_back(model);
  }
  return true;
}

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

// VHE DATA READER

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

class VHEData {
public:
  VHEData(const std::string &filename, double deltae_over_e = 0.0);
  unsigned numData() const { return fEnergy.size(); }
  const double *e() const { return &(*fEnergy.begin()); }
  const double *e2dNdE() const { return &(*fE2dNdE.begin()); }
  const double *deltaE2dNdE() const { return &(*fDeltaE2dNdE.begin()); }
  bool redshift(double &z) const;
  std::string name() const;
  bool parameter(const std::string &key, std::string &val) const;

private:
  double fEMultiplier;
  std::string fName;
  std::map<std::string, std::string> fParameters;
  std::vector<double> fEnergy;
  std::vector<double> fE2dNdE;
  std::vector<double> fDeltaE2dNdE;
  bool fCommentToParameters;
  bool getline(std::istream &stream, std::string &line);
};

#endif // EBLMODELSFILE_HPP
