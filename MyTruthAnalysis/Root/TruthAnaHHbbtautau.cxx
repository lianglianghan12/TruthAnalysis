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

  std::cout << "-----------------------> it is here" << std::endl;
  std::cout << "in execute, runNumber = " << eventInfo->runNumber() << ", eventNumber = " << eventInfo->eventNumber() << std::endl;

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

  std::cout << "m_fMCWeight is " << m_fMCWeight << std::endl;

  // retrieve the TruthParticles container.
  const xAOD::TruthParticleContainer *truthParticleContainer = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthParticleContainer, "TruthParticles"));
  vector<const xAOD::TruthParticle *> finalStateParticleVec{};
  ANA_MSG_DEBUG("the size of the truthParticleContainer is: "<<truthParticleContainer->size());
  std::cout << "the size of the truthParticleContainer is: "<<truthParticleContainer->size() << std::endl;
  for (std::size_t i = 0; i < truthParticleContainer->size(); i++)
  {
    if (truthParticleContainer->at(i)->status() == 1)
    {
      finalStateParticleVec.push_back(truthParticleContainer->at(i));
      //std::cout << "the PDG ID of the each particle in the truthParticle container: " << truthParticleContainer->at(i)->pdgId() << std::endl;
      //std::cout << "the eta of the particle is: " << truthParticleContainer->at(i)->eta() << std::endl;
      //std::cout << "the number of children of the particle is: " << truthParticleContainer->at(i)->nChildren() << std::endl;
      ANA_MSG_DEBUG("the PDG ID of the each particle in the truthParticle container: " << truthParticleContainer->at(i)->pdgId());
      ANA_MSG_DEBUG("the eta of the particle is: " << truthParticleContainer->at(i)->eta());
      ANA_MSG_DEBUG("the number of children of the particle is: " << truthParticleContainer->at(i)->nChildren());
    }
  }
//l  if (incomingQuarkVec.size() != 2) return StatusCode::FAILURE; // actually, I found that only 2 incoming partons are stored for all events, this statement is used to avoid exceptional events.
//l  if (incomingQuarkVec[0]->nChildren() != 4) return StatusCode::SUCCESS;
//l
//l  // retrieve jet container where the jets are "reconstructed" with the "AntiKt4TruthDressedWZJets" algorithm
//l  const xAOD::JetContainer *jets = nullptr;
//l  ANA_CHECK(evtStore()->retrieve(jets, "AntiKt4TruthDressedWZJets"));  
//l  APPLYCUT(true, "Initial");
//l//  APPLYCUT(jets->size() > 5, "Number of truth jets")
//l  m_nJets = jets->size();
//l  
//l  // Vector and map definition for housing objects. Specifically, the vbfTagJetVec is for VBF-tag jets candidates while the BjetsJetVec is for b-jets
//l  std::map<std::string, const xAOD::TruthParticle *> higgsMap{};
//l//  vector<const xAOD::Jet *> vbfTagJetVec{};
//l//  vector<const xAOD::Jet *> BjetsJetVec{};
//l  vector<const xAOD::TruthParticle *> outgoingQuark0Vec{};
//l  vector<const xAOD::TruthParticle *> outgoingQuark1Vec{};
//l
//l  // fill up the outgoing quark vectors with the quarks separately from the eligible children of the 2 incoming quarks. Idea: put the eligible children of the 2 incoming quarks separately into 2 vectors.
//l  outgoingQuark0Vec.push_back( getFinal( incomingQuarkVec[0]->child(2) ) );
//l  outgoingQuark1Vec.push_back( getFinal( incomingQuarkVec[0]->child(3) ) );
//l  ANA_MSG_DEBUG("the vbf quark0 is: " << outgoingQuark0Vec[0]->pdgId());
//l  ANA_MSG_DEBUG("the vbf quark1 is: " << outgoingQuark1Vec[0]->pdgId());
//l
//l  // fill up the higgsMap map with the higgs from the truthEvent container. Idea: loop over all the particles in the truthEvent container, put the particle with the pdgId of 25 and only 2 tau children into higgsMap as well as the particle with the pdgId of 25 and only 2 b children
//l  for (std::size_t i = 0; i < 2; i++)
//l  {
//l    const xAOD::TruthParticle *particle = getFinal( incomingQuarkVec[0]->child(i) );
//l
//l    // fetch Higgs
//l    vector<unsigned> trash_tautau_idx;
//l    if (particle->pdgId() == 25 && particle->nChildren() == 2 && hasChild(particle, 15, trash_tautau_idx))
//l    {
//l      higgsMap.insert({"Htautau", particle});
//l    }
//l
//l    vector<unsigned> trash_bb_idx;
//l    if (particle->pdgId() == 25 && particle->nChildren() == 2 && hasChild(particle, 5, trash_tautau_idx))
//l    {
//l      higgsMap.insert({"Hbb", particle});
//l    }
//l  }
//l  if (higgsMap.size() != 2) return StatusCode::FAILURE;
//l
//l  // particles
//l  const xAOD::TruthParticle *tau0 = getFinal(higgsMap["Htautau"]->child(0));
//l  const xAOD::TruthParticle *tau1 = getFinal(higgsMap["Htautau"]->child(1));
//l  const xAOD::TruthParticle* b0   = getFinal(higgsMap["Hbb"]->child(0));
//l  const xAOD::TruthParticle* b1   = getFinal(higgsMap["Hbb"]->child(1));
//l  // const xAOD::TruthParticle *b0 = higgsMap["Hbb"]->child(0);
//l  // const xAOD::TruthParticle *b1 = higgsMap["Hbb"]->child(1);
//l
//l  // kinematics
//l  TLorentzVector tau0_p4 = tau0->p4(), tau1_p4 = tau1->p4();
//l  TLorentzVector tauvis0_p4 = tauVisP4(tau0), tauvis1_p4 = tauVisP4(tau1);
//l  TLorentzVector b0_p4 = b0->p4(), b1_p4 = b1->p4();
//l
//l  ANA_MSG_DEBUG("Htautau : " << higgsMap["Htautau"]->child(0)->pdgId() << ", " << higgsMap["Htautau"]->child(1)->pdgId());
//l  ANA_MSG_DEBUG("Hbb     : " << higgsMap["Hbb"]->child(0)->pdgId() << ", " << higgsMap["Hbb"]->child(1)->pdgId());
//l
//l  APPLYCUT(isGoodEvent(), "Empty cut for testing");
//l  APPLYCUT(isOS(tau0, tau1) && isOS(b0, b1), "OS Charge");
//l  APPLYCUT(isGoodTau(tau0, 20., 2.5) && isGoodTau(tau1, 20., 2.5), "Tau Preselection");//Tau candidate reco requirement
//l  APPLYCUT(isGoodB(b0, 20., 2.4) && isGoodB(b1, 20., 2.4), "B-jet preselection");//B-jet candidate reco requirement
//l  APPLYCUT(isNotOverlap(b0, b1, tau0, tau1, 0.2), "b-tau overlap removal");
//l
//l  // mimic single tau trigger selection
//l  bool STT = isGoodTau(tau0, 100., 2.5) && isGoodB(b0, 45., 2.4);
//l
//l  // mimic di-tau trigger selection
//l  bool DTT = isGoodTau(tau0, 40., 2.5) && isGoodTau(tau1, 30., 2.5) && isGoodB(b0, 80., 2.4);
//l
//l  // event must pass single tau trigger or di-tau trigger
//l  APPLYCUT(STT || DTT, "Trigger selection (TO CHECK)");
//l  APPLYCUT((tau0_p4 + tau1_p4).M() > 60 * GeV, "Di-tau mass selection");//to reduce Drell-Yan in hadhad channel
//l
//l  // 4-momenta
//l  m_fTau0_pt = tau0_p4.Pt() / GeV;
//l  m_fTau0_phi = tau0_p4.Phi();
//l  m_fTau0_eta = tau0_p4.Eta();
//l  m_fTau1_pt = tau1_p4.Pt() / GeV;
//l  m_fTau1_phi = tau1_p4.Phi();
//l  m_fTau1_eta = tau1_p4.Eta();
//l
//l  m_fTauVis0_pt = tauvis0_p4.Pt() / GeV;
//l  m_fTauVis0_phi = tauvis0_p4.Phi();
//l  m_fTauVis0_eta = tauvis0_p4.Eta();
//l  m_fTauVis1_pt = tauvis1_p4.Pt() / GeV;
//l  m_fTauVis1_phi = tauvis1_p4.Phi();
//l  m_fTauVis1_eta = tauvis1_p4.Eta();
//l
//l  m_fB0_pt = b0_p4.Pt() / GeV;
//l  m_fB0_phi = b0_p4.Phi();
//l  m_fB0_eta = b0_p4.Eta();
//l  m_fB1_pt = b1_p4.Pt() / GeV;
//l  m_fB1_phi = b1_p4.Phi();
//l  m_fB1_eta = b1_p4.Eta();
//l
//l  // deltaRs
//l  m_fDeltaR_TauTau = tau0_p4.DeltaR(tau1_p4);
//l  m_fDeltaR_TauVisTauVis = tauvis0_p4.DeltaR(tauvis1_p4);
//l  m_fDeltaR_BB = b0_p4.DeltaR(b1_p4);
//l  m_fDeltaR_BB_TauTau = (b0_p4 + b1_p4).DeltaR(tau0_p4 + tau1_p4);
//l
//l  // Higgs Pt
//l  m_fPtBB = (b0_p4 + b1_p4).Pt() / GeV;
//l  m_fPtTauTau = (tau0_p4 + tau1_p4).Pt() / GeV;
//l
//l  // invariant masses
//l  m_fMTauTau = (tau0_p4 + tau1_p4).M() / GeV;
//l  m_fMTauVisTauVis = (tauvis0_p4 + tauvis1_p4).M() / GeV;
//l  m_fMBB = (b0_p4 + b1_p4).M() / GeV;
//l  m_fMHH = (tau0_p4 + tau1_p4 + b0_p4 + b1_p4).M() / GeV;
//l
//l  // select the vbf truth quark
//l  if (outgoingQuark0Vec.size() != 1 || outgoingQuark1Vec.size() != 1) 
//l  {
//l    return StatusCode::FAILURE;
//l  }
//l  else
//l  {
//l    const xAOD::TruthParticle *VTjet0 = outgoingQuark0Vec[0];
//l    const xAOD::TruthParticle *VTjet1 = outgoingQuark1Vec[0];
//l    APPLYCUT( (VTjet0->eta())*(VTjet1->eta()) < 0, "eta OS cut");
//l    ANA_MSG_DEBUG("--------------------------------------Pt is " << VTjet0->p4().Pt());
//l    ANA_MSG_DEBUG("--------------------------------------Eta is " << VTjet0->p4().Eta());
//l    ANA_MSG_DEBUG("--------------------------------------Phi is " << VTjet0->p4().Phi());
//l
//l    TLorentzVector VTjet0_p4, VTjet1_p4;
//l    VTjet0_p4 = VTjet0->p4();
//l    VTjet1_p4 = VTjet1->p4();
//l    m_fVTjet0_pt = VTjet0_p4.Pt() / GeV;
//l    m_fVTjet0_phi = VTjet0_p4.Phi();
//l    m_fVTjet0_eta = VTjet0_p4.Eta();
//l    m_fVTjet1_pt = VTjet1_p4.Pt() / GeV;
//l    m_fVTjet1_phi = VTjet1_p4.Phi();
//l    m_fVTjet1_eta = VTjet1_p4.Eta();
//l    m_fDeltaR_VTjetVTjet = VTjet0_p4.DeltaR(VTjet1_p4);
//l    m_fDeltaEta_VTjetVTjet = abs(m_fVTjet0_eta - m_fVTjet1_eta);
//l    m_fMVTjetVTjet = (VTjet0_p4 + VTjet1_p4).M() / GeV;
//l  }
//l
//l  // apply addition cuts for VBF tag jets
//l//l  APPLYCUT(m_fMVTjetVTjet > 1000, "Mjj cut");
//l  APPLYCUT(m_fDeltaEta_VTjetVTjet > 3, "DeltaEta_jj");
//l//l  APPLYCUT(m_fVTjet0_eta * m_fVTjet1_eta < 0, "Eta*EtaOS_jj");

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
  m_cTree->Branch("NJets", &m_nJets);
  m_cTree->Branch("NFatJets", &m_nFatJets);
  m_cTree->Branch("Tau0_pt", &m_fTau0_pt);
  m_cTree->Branch("Tau1_pt", &m_fTau1_pt);
  m_cTree->Branch("Tau0_phi", &m_fTau0_phi);
  m_cTree->Branch("Tau1_phi", &m_fTau1_phi);
  m_cTree->Branch("Tau0_eta", &m_fTau0_eta);
  m_cTree->Branch("Tau1_eta", &m_fTau1_eta);
  m_cTree->Branch("TauVis0_pt", &m_fTauVis0_pt);
  m_cTree->Branch("TauVis1_pt", &m_fTauVis1_pt);
  m_cTree->Branch("TauVis0_phi", &m_fTauVis0_phi);
  m_cTree->Branch("TauVis1_phi", &m_fTauVis1_phi);
  m_cTree->Branch("TauVis0_eta", &m_fTauVis0_eta);
  m_cTree->Branch("TauVis1_eta", &m_fTauVis1_eta);
  m_cTree->Branch("B0_pt", &m_fB0_pt);
  m_cTree->Branch("B1_pt", &m_fB1_pt);
  m_cTree->Branch("B0_phi", &m_fB0_phi);
  m_cTree->Branch("B1_phi", &m_fB1_phi);
  m_cTree->Branch("B0_eta", &m_fB0_eta);
  m_cTree->Branch("B1_eta", &m_fB1_eta);
  m_cTree->Branch("Bjet0_pt", &m_fBjet0_pt);
  m_cTree->Branch("Bjet1_pt", &m_fBjet1_pt);
  m_cTree->Branch("Bjet0_phi", &m_fBjet0_phi);
  m_cTree->Branch("Bjet1_phi", &m_fBjet1_phi);
  m_cTree->Branch("Bjet0_eta", &m_fBjet0_eta);
  m_cTree->Branch("Bjet1_eta", &m_fBjet1_eta);
  m_cTree->Branch("DiBjet_pt", &m_fDiBjet_pt);
  m_cTree->Branch("DiBjet_m", &m_fDiBjet_m);
  m_cTree->Branch("DiBjet_phi", &m_fDiBjet_phi);
  m_cTree->Branch("DiBjet_eta", &m_fDiBjet_eta);
  m_cTree->Branch("DeltaR_BB", &m_fDeltaR_BB);
  m_cTree->Branch("DeltaR_BjetBjet", &m_fDeltaR_BjetBjet);
  m_cTree->Branch("DeltaR_TauTau", &m_fDeltaR_TauTau);
  m_cTree->Branch("DeltaR_TauVisTauVis", &m_fDeltaR_TauVisTauVis);
  m_cTree->Branch("DeltaR_BB_TauTau", &m_fDeltaR_BB_TauTau);
  m_cTree->Branch("DeltaR_BjetBjet_TauVisTauVis", &m_fDeltaR_BjetBjet_TauVisTauVis);
  m_cTree->Branch("DeltaR_DiBjet_TauVisTauVis", &m_fDeltaR_DiBjet_TauVisTauVis);
  m_cTree->Branch("MBB", &m_fMBB);
  m_cTree->Branch("MTauTau", &m_fMTauTau);
  m_cTree->Branch("MBjetBjet", &m_fMBjetBjet);
  m_cTree->Branch("MTauVisTauVis", &m_fMTauVisTauVis);
  m_cTree->Branch("MHH", &m_fMHH);
  m_cTree->Branch("PtBB", &m_fPtBB);
  m_cTree->Branch("PtTauTau", &m_fPtTauTau);
  m_cTree->Branch("MCWeight", &m_fMCWeight);
  m_cTree->Branch("Channel", &m_nChannel);
  m_cTree->Branch("VTjet0_pt", &m_fVTjet0_pt);
  m_cTree->Branch("VTjet0_phi", &m_fVTjet0_phi);
  m_cTree->Branch("VTjet0_eta", &m_fVTjet0_eta);
  m_cTree->Branch("VTjet1_pt", &m_fVTjet1_pt);
  m_cTree->Branch("VTjet1_phi", &m_fVTjet1_phi);
  m_cTree->Branch("VTjet1_eta", &m_fVTjet1_eta);
  m_cTree->Branch("DeltaR_VTjetVTjet", &m_fDeltaR_VTjetVTjet);
  m_cTree->Branch("DeltaEta_VTjetVTjet", &m_fDeltaEta_VTjetVTjet);
  m_cTree->Branch("MVTjetVTjet", &m_fMVTjetVTjet);
}

void TruthAnaHHbbtautau::resetBranches()
{
  m_nEventNumber = 0; // DONE
  m_nRunNumber = 0; // DONE
  m_nJets = 0; // DONE
  m_nFatJets = 0; // DONE
  m_fTau0_pt = 0; // DONE
  m_fTau1_pt = 0; // DONE
  m_fTau0_phi = 0; // DONE
  m_fTau1_phi = 0; // DONE
  m_fTau0_eta = 0; // DONE
  m_fTau1_eta = 0; // DONE
  m_fTauVis0_pt = 0; // DONE
  m_fTauVis1_pt = 0; // DONE
  m_fTauVis0_phi = 0; // DONE
  m_fTauVis1_phi = 0; // DONE
  m_fTauVis0_eta = 0; // DONE
  m_fTauVis1_eta = 0; // DONE
  m_fB0_pt = 0; // DONE
  m_fB1_pt = 0; // DONE
  m_fB0_phi = 0; // DONE
  m_fB1_phi = 0; // DONE
  m_fB0_eta = 0; // DONE
  m_fB1_eta = 0; // DONE
  m_fBjet0_pt = 0; // DONE
  m_fBjet1_pt = 0; // DONE
  m_fBjet0_phi = 0; // DONE
  m_fBjet1_phi = 0; // DONE
  m_fBjet0_eta = 0; // DONE
  m_fBjet1_eta = 0; // DONE
  m_fDiBjet_pt = 0; // DONE
  m_fDiBjet_m = 0; // DONE
  m_fDiBjet_phi = 0; // DONE
  m_fDiBjet_eta = 0; // DONE
  m_fDeltaR_BB = 0; // DONE
  m_fDeltaR_BjetBjet = 0; // DONE
  m_fDeltaR_TauTau = 0; // DONE
  m_fDeltaR_TauVisTauVis = 0; // DONE
  m_fDeltaR_BB_TauTau = 0; // DONE
  m_fDeltaR_BjetBjet_TauVisTauVis = 0; // DONE
  m_fDeltaR_DiBjet_TauVisTauVis = 0; // DONE
  m_fMBB = 0; // DONE
  m_fMTauTau = 0; // DONE
  m_fMBjetBjet = 0; // DONE
  m_fMTauVisTauVis = 0; // DONE
  m_fMHH = 0; // DONE
  m_fPtBB = 0; // DONE
  m_fPtTauTau = 0; // DONE
  m_fMCWeight = 0; // DONE
  m_nChannel = 0; // DONE
  m_fVTjet0_pt = 0;                // DONE
  m_fVTjet0_phi = 0;               // DONE
  m_fVTjet0_eta = 0;               // DONE
  m_fVTjet1_pt = 0;                // DONE
  m_fVTjet1_phi = 0;               // DONE
  m_fVTjet1_eta = 0;               // DONE
  m_fDeltaR_VTjetVTjet = 0;        // DONE
  m_fDeltaEta_VTjetVTjet = 0;      // DONE
  m_fMVTjetVTjet = 0;              // DONE
}
