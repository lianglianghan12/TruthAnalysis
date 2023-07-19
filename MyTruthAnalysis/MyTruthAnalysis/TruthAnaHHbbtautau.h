#ifndef MyTruthAnalysis_TruthAnaHHbbtautau_H
#define MyTruthAnalysis_TruthAnaHHbbtautau_H

// Base class
#include <AnaAlgorithm/AnaAlgorithm.h>

// My class
#include "MyTruthAnalysis/Cutflow.h"
#include "MyTruthAnalysis/TruthAnaBase.h"

// std
#include <memory>

class TTree;

enum class CHAN { UNKNOWN=0, RESOLVED, BOOSTED }; //!

class TruthAnaHHbbtautau : public TruthAnaBase
{
public:
  // this is a standard algorithm constructor
  TruthAnaHHbbtautau(const std::string &name, ISvcLocator *pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  void initBranches();
  void resetBranches();

private:
  TTree *m_cTree = nullptr;          //!
  CHAN m_eChannel = CHAN::UNKNOWN;   //!

private:
  unsigned long long m_nEventNumber; //!
  unsigned long long m_nRunNumber;   //!

  double m_fMCWeight;                //!

  unsigned long long m_nChannel;     //!
};

#endif
