// AsgTools
#include <AsgTools/MessageCheck.h>
#include <AsgTools/MsgLevel.h>

// xAOD
#include <xAODEventInfo/EventInfo.h>
#include <xAODTruth/TruthEvent.h>
#include <xAODTruth/TruthEventContainer.h>
#include <xAODTruth/TruthParticle.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthVertex.h>
#include <xAODJet/JetContainer.h>
#include <xAODJet/Jet.h>

// ROOT
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1.h>


// My headers
#include "MyTruthAnalysis/TruthAnaHHbbtautau.h"
#include "MyTruthAnalysis/HelperFunctions.h"

// std
#include <map>
#include <memory>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cassert>

using namespace TruthAna;

using std::vector;

TruthAnaHHbbtautau::TruthAnaHHbbtautau(const std::string &name,
                                       ISvcLocator *pSvcLocator)
    : TruthAnaBase(name, pSvcLocator)
{
}

StatusCode TruthAnaHHbbtautau::initialize()
{
  ANA_MSG_INFO("Initializing ...");
  ANA_CHECK( book( TTree("MyTree", "truth analysis tree") ) );
  m_cTree = tree("MyTree");
  m_cTree->SetMaxTreeSize(500'000'000);
  initBranches();

  return StatusCode::SUCCESS;
}

StatusCode TruthAnaHHbbtautau::execute()
{
  resetBranches();

  // retrieve the eventInfo object from the event store
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK(evtStore()->retrieve(eventInfo, "EventInfo"));
  ANA_MSG_DEBUG("in execute, runNumber = " << eventInfo->runNumber() << ", eventNumber = " << eventInfo->eventNumber());
  m_nRunNumber = eventInfo->runNumber();
  m_nEventNumber = eventInfo->eventNumber();
  m_nChannel = static_cast<unsigned long long>(m_eChannel);

  // retrieve truthEvent Container
  const xAOD::TruthEventContainer *truthEventContainer = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthEventContainer, "TruthEvents"));
  const xAOD::TruthEvent *truthEvent = nullptr;
  if (truthEventContainer->size() == 1)
    truthEvent = truthEventContainer->at(0);
  else
  {
    ANA_MSG_WARNING("in execute, no truth event container!");
    return StatusCode::SUCCESS;
  }
  const vector<float> weights = truthEvent->weights();
  m_fMCWeight = weights[0];

  // retrieve the TruthParticles container.
  const xAOD::TruthParticleContainer *truthParticleContainer = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthParticleContainer, "TruthParticles"));
  vector<const xAOD::TruthParticle *> finalStateParticleVec{};
  ANA_MSG_DEBUG("the size of the truthParticleContainer is: "<<truthParticleContainer->size());
  for (std::size_t i = 0; i < truthParticleContainer->size(); i++)
  {
    if (truthParticleContainer->at(i)->status() == 1)
    {
      finalStateParticleVec.push_back(truthParticleContainer->at(i));
      ANA_MSG_DEBUG("the PDG ID of the each particle in the truthParticle container: " << truthParticleContainer->at(i)->pdgId());
      ANA_MSG_DEBUG("the eta of the particle is: " << truthParticleContainer->at(i)->eta());
      ANA_MSG_DEBUG("the number of children of the particle is: " << truthParticleContainer->at(i)->nChildren());
    }
  }

  m_cTree->Fill();

  return StatusCode::SUCCESS;
}

StatusCode TruthAnaHHbbtautau::finalize()
{
  ANA_MSG_INFO("Finalizing ...");
  m_cCutflow->print();
  return StatusCode::SUCCESS;
}

// ----------------------------------------------------------------------------
// Function for trees
// ----------------------------------------------------------------------------

void TruthAnaHHbbtautau::initBranches()
{
  m_cTree->Branch("EventNumber", &m_nEventNumber);
  m_cTree->Branch("RunNumber", &m_nRunNumber);
  m_cTree->Branch("MCWeight", &m_fMCWeight);
  m_cTree->Branch("Channel", &m_nChannel);
}

void TruthAnaHHbbtautau::resetBranches()
{
  m_nEventNumber = 0; // DONE
  m_nRunNumber = 0; // DONE
  m_fMCWeight = 0; // DONE
  m_nChannel = 0; // DONE
}
