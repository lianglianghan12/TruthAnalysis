// AsgTools
#include <AsgTools/MessageCheck.h>

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
  vector<const xAOD::TruthParticle *> incomingQuarkVec{};
  ANA_MSG_DEBUG("the size of the truthParticleContainer is: "<<truthParticleContainer->size());
  for (std::size_t i = 0; i < truthParticleContainer->size(); i++)
  {
    if (truthParticleContainer->at(i)->status() == 4)
    {
      incomingQuarkVec.push_back(truthParticleContainer->at(i));
      ANA_MSG_DEBUG("the PDG ID of the each particle in the truthParticle container: " << truthParticleContainer->at(i)->pdgId());
    }
  }
  if (incomingQuarkVec.size() != 2) return StatusCode::FAILURE; // actually, I found that only 2 incoming partons are stored for all events, this statement is used to avoid exceptional events.

  // retrieve jet container where the jets are "reconstructed" with the "AntiKt4TruthDressedWZJets" algorithm
  const xAOD::JetContainer *jets = nullptr;
  ANA_CHECK(evtStore()->retrieve(jets, "AntiKt4TruthDressedWZJets"));  
  APPLYCUT(true, "Initial");
//  APPLYCUT(jets->size() > 5, "Number of truth jets")
  m_nJets = jets->size();
  
  // Vector and map definition for housing objects. Specifically, the vbfTagJetVec is for VBF-tag jets candidates while the BjetsJetVec is for b-jets
  std::map<std::string, const xAOD::TruthParticle *> higgsMap{};
  vector<const xAOD::Jet *> vbfTagJetVec{};
  vector<const xAOD::Jet *> BjetsJetVec{};
  vector<const xAOD::TruthParticle *> outgoingQuark0Vec{};
  vector<const xAOD::TruthParticle *> outgoingQuark1Vec{};

  // fill up the outgoing quark vectors with the quarks separately from the eligible children of the 2 incoming quarks. Idea: put the eligible children of the 2 incoming quarks separately into 2 vectors.
  for (std::size_t iChild = 0; iChild < incomingQuarkVec[0]->nChildren(); iChild++)
  {
    const xAOD::TruthParticle *childG = getFinal(incomingQuarkVec[0]->child(iChild));
//    ANA_MSG_DEBUG("The number of children for the 1st incoming parton: " << incomingQuarkVec[0]->nChildren());
//    ANA_MSG_DEBUG("The Eta of the " << iChild << "th particle is " << childG->p4().Eta());
//    ANA_MSG_DEBUG("The Phi of the " << iChild << "th particle is " << childG->p4().Phi()); 
//    ANA_MSG_DEBUG("The Pt of the " << iChild << "th particle is " << (childG->p4().Pt()) / GeV);
//    ANA_MSG_DEBUG("The M of the " << iChild << "th particle is " << (childG->m()) / GeV);
//    ANA_MSG_DEBUG("The pdgId of the " << iChild << "th particle is " << childG->pdgId());
//    if (childG->absPdgId() > -5 && childG->absPdgId() < 5 && childG->absPdgId() != 0)
    if (childG->pdgId() == 82)
    {
      outgoingQuark0Vec.push_back(childG);
//      ANA_MSG_DEBUG("The number of children for the 1st incoming parton: " << incomingQuarkVec[0]->nChildren());
//      ANA_MSG_DEBUG("The Eta of the " << iChild << "th particle is " << childG->p4().Eta()); 
//      ANA_MSG_DEBUG("The Phi of the " << iChild << "th particle is " << childG->p4().Phi()); 
//      ANA_MSG_DEBUG("The Pt of the " << iChild << "th particle is " << (childG->p4().Pt()) / GeV);
//      ANA_MSG_DEBUG("The M of the " << iChild << "th particle is " << (childG->m()) / GeV);
//      ANA_MSG_DEBUG("The pdgId of the " << iChild << "th particle is " << childG->pdgId());
    }
  }
  for (std::size_t iChild = 0; iChild < incomingQuarkVec[1]->nChildren(); iChild++)
  {
    const xAOD::TruthParticle *childH = getFinal(incomingQuarkVec[1]->child(iChild));
//    ANA_MSG_DEBUG("The number of children for the 2nd incoming parton: " << incomingQuarkVec[1]->nChildren());
//    ANA_MSG_DEBUG("The Eta of the " << iChild << "th particle is " << childH->p4().Eta());
//    ANA_MSG_DEBUG("The Phi of the " << iChild << "th particle is " << childH->p4().Phi());
//    ANA_MSG_DEBUG("The Pt of the " << iChild << "th particle is " << (childH->p4().Pt()) / GeV);
//    ANA_MSG_DEBUG("The M of the " << iChild << "th particle is " << (childH->m()) / GeV);
//    ANA_MSG_DEBUG("The pdgId of the " << iChild << "th particle is " << childH->pdgId());
//    if (childH->absPdgId() > -5 && childH->absPdgId() < 5 && childH->absPdgId() != 0)
    if (childH->pdgId() == 82)
    {
      outgoingQuark1Vec.push_back(childH);
//      ANA_MSG_DEBUG("The number of children for the 2nd incoming parton: " << incomingQuarkVec[1]->nChildren());
//      ANA_MSG_DEBUG("The Eta of the " << iChild << "th particle is " << childH->p4().Eta()); 
//      ANA_MSG_DEBUG("The Phi of the " << iChild << "th particle is " << childH->p4().Phi()); 
//      ANA_MSG_DEBUG("The Pt of the " << iChild << "th particle is " << (childH->p4().Pt()) / GeV);
//      ANA_MSG_DEBUG("The M of the " << iChild << "th particle is " << (childH->m()) / GeV);
//      ANA_MSG_DEBUG("The pdgId of the " << iChild << "th particle is " << childH->pdgId());
    }
  }
//  if (outgoingQuark0Vec.size() > 0 && outgoingQuark1Vec.size() > 0)
//  {
//    ANA_MSG_DEBUG("the Mjj 0 0 summation is:" << ((outgoingQuark0Vec[0]->p4() + outgoingQuark1Vec[0]->p4()).M()) / GeV);
//  }
  if (outgoingQuark0Vec.size() != 1 || outgoingQuark1Vec.size() != 1)
  {
    return StatusCode::FAILURE;
  }
  ANA_MSG_DEBUG("The number of the child of the outgoing parton0: "<< outgoingQuark0Vec[0]->nChildren());
  for (std::size_t iChild = 0; iChild < outgoingQuark0Vec[0]->nChildren(); iChild++)
  {
    ANA_MSG_DEBUG("the pdgid of the " << iChild << " th child is: " << getFinal(outgoingQuark0Vec[0]->child(iChild))->pdgId());
  }
  ANA_MSG_DEBUG("The number of the child of the outgoing parton1: "<< outgoingQuark1Vec[0]->nChildren());
  for (std::size_t iChild = 0; iChild < outgoingQuark1Vec[0]->nChildren(); iChild++)
  {
    ANA_MSG_DEBUG("the pdgid of the " << iChild << " th child is: " << getFinal(outgoingQuark1Vec[0]->child(iChild))->pdgId());
  }

  // fill up the higgsMap map with the higgs from the truthEvent container. Idea: loop over all the particles in the truthEvent container, put the particle with the pdgId of 25 and only 2 tau children into higgsMap as well as the particle with the pdgId of 25 and only 2 b children
  bool found_htautau = false;
  bool found_hbb = false;
  ANA_MSG_DEBUG("The number of the particles in the truth event : " << truthEvent->nTruthParticles());
  for (std::size_t i = 0; i < truthEvent->nTruthParticles(); i++)
  {
    const xAOD::TruthParticle *particle = truthEvent->truthParticle(i);

    // Some debug output
//    if (i == 0)
//    {
//      ANA_MSG_DEBUG("Particle info: ");
//      ANA_MSG_DEBUG(" - Barcode: " << particle->barcode());
//      ANA_MSG_DEBUG(" - PDG ID : " << particle->pdgId());
//      ANA_MSG_DEBUG(" - Pt     : " << particle->pt());
//    }

    // fetch Higgs
    vector<unsigned> trash_tautau_idx;
    if (!found_htautau && particle->pdgId() == 25 && particle->nChildren() == 2 && hasChild(particle, 15, trash_tautau_idx))
    {
//      found_htautau = true;
      higgsMap.insert({"Htautau", particle});
    }

    vector<unsigned> trash_bb_idx;
    if (!found_hbb && particle->pdgId() == 25 && particle->nChildren() == 2 && hasChild(particle, 5, trash_tautau_idx))
    {
//      found_hbb = true;
      higgsMap.insert({"Hbb", particle});
    }

//    if (found_hbb && found_htautau) break;
  }
  if (higgsMap.size() != 2) return StatusCode::FAILURE;

  // fill up the vbfTagJetVec vector with the vbf-tag jets from the jets container. Idea: except for the more than 1 tau jets and more than 1 b-jets, select all the remaining jets as the vbf-jet candidates.
//  ANA_MSG_DEBUG("The number of jets in the jets container : " << jets->size());
//  vector<int> vbfJetTag_idx{};
//  vector<int> JetsLableID{};
//  for (std::size_t i = 0; i < jets->size(); i++)
//  {
//    ANA_MSG_DEBUG("Jet truth flavour info: ");
//    ANA_MSG_DEBUG(" - PartonTruthLabelID = " << jets->at(i)->auxdata<int>("PartonTruthLabelID"));
//    ANA_MSG_DEBUG(" - HadronConeExclTruthLabelID = " << jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID"));
//    ANA_MSG_DEBUG(" - ConeTruthLabelID = " << jets->at(i)->auxdata<int>("ConeTruthLabelID"));
//    JetsLableID.push_back(jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID"));
//  }
//  if ( std::count(JetsLableID.begin(), JetsLableID.end(), 15) > 1 && std::count(JetsLableID.begin(), JetsLableID.end(), 5) > 1 )
//  {
//    for (std::size_t i = 0; i < jets->size(); i++)
//    {
//      if (jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID") == 15) //Check if it is TauJet, if so, pass
//      {
//        continue;
//      }
//      else if (jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID") == 5) //Check if it is b-jet, if so, pass
//      {
//        continue;
//      }
//      else
//      {
//        vbfTagJetVec.push_back(jets->at(i));
//        vbfJetTag_idx.push_back(i);
//      }
//    }
//  }
//  else
//  {
//    return StatusCode::SUCCESS;   
//  }
//  ANA_MSG_DEBUG("vbfTagJetVec size is: " << vbfTagJetVec.size());

  // fill up the BjetsJetVec vector with the b-jets from the jets container. Idea: just take the jets with the pdgid of 5 as the b-jets
//  ANA_MSG_DEBUG("The number of jets in the jets container : " << jets->size());
//  vector<int> BJetTag_idx{};
//  for (std::size_t i = 0; i < jets->size(); i++)
//  {
//    ANA_MSG_DEBUG("Jet truth flavour info: ");
//    ANA_MSG_DEBUG(" - PartonTruthLabelID = " << jets->at(i)->auxdata<int>("PartonTruthLabelID"));
//    ANA_MSG_DEBUG(" - HadronConeExclTruthLabelID = " << jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID"));
//    ANA_MSG_DEBUG(" - ConeTruthLabelID = " << jets->at(i)->auxdata<int>("ConeTruthLabelID"));
//    if (jets->at(i)->auxdata<int>("HadronConeExclTruthLabelID") == 5)
//    { // isBJet -> TruthFlavor == 5
//      BjetsJetVec.push_back(jets->at(i));
//      BJetTag_idx.push_back(i);
//    }
//  }

  // particles
  const xAOD::TruthParticle *tau0 = getFinal(higgsMap["Htautau"]->child(0));
  const xAOD::TruthParticle *tau1 = getFinal(higgsMap["Htautau"]->child(1));
  const xAOD::TruthParticle* b0   = getFinal(higgsMap["Hbb"]->child(0));
  const xAOD::TruthParticle* b1   = getFinal(higgsMap["Hbb"]->child(1));
  // const xAOD::TruthParticle *b0 = higgsMap["Hbb"]->child(0);
  // const xAOD::TruthParticle *b1 = higgsMap["Hbb"]->child(1);

  // kinematics
  TLorentzVector tau0_p4 = tau0->p4(), tau1_p4 = tau1->p4();
  TLorentzVector tauvis0_p4 = tauVisP4(tau0), tauvis1_p4 = tauVisP4(tau1);
  TLorentzVector b0_p4 = b0->p4(), b1_p4 = b1->p4();

  ANA_MSG_DEBUG("Htautau : " << higgsMap["Htautau"]->child(0)->pdgId() << ", " << higgsMap["Htautau"]->child(1)->pdgId());
  ANA_MSG_DEBUG("Hbb     : " << higgsMap["Hbb"]->child(0)->pdgId() << ", " << higgsMap["Hbb"]->child(1)->pdgId());

//  APPLYCUT(isGoodEvent(), "Empty cut for testing");
//  APPLYCUT(isOS(tau0, tau1) && isOS(b0, b1), "OS Charge");
//  APPLYCUT(isGoodTau(tau0, 20., 2.5) && isGoodTau(tau1, 20., 2.5), "Tau Preselection");//Tau candidate reco requirement
//  APPLYCUT(isGoodB(b0, 20., 2.4) && isGoodB(b1, 20., 2.4), "B-jet preselection");//B-jet candidate reco requirement
//  APPLYCUT(isNotOverlap(b0, b1, tau0, tau1, 0.2), "b-tau overlap removal");

  // mimic single tau trigger selection
  bool STT = isGoodTau(tau0, 100., 2.5) && isGoodB(b0, 45., 2.4);

  // mimic di-tau trigger selection
  bool DTT = isGoodTau(tau0, 40., 2.5) && isGoodTau(tau1, 30., 2.5) && isGoodB(b0, 80., 2.4);

  // event must pass single tau trigger or di-tau trigger
//  APPLYCUT(STT || DTT, "Trigger selection (TO CHECK)");
//  APPLYCUT((tau0_p4 + tau1_p4).M() > 60 * GeV, "Di-tau mass selection");//to reduce Drell-Yan in hadhad channel

  // 4-momenta
  m_fTau0_pt = tau0_p4.Pt() / GeV;
  m_fTau0_phi = tau0_p4.Phi();
  m_fTau0_eta = tau0_p4.Eta();
  m_fTau1_pt = tau1_p4.Pt() / GeV;
  m_fTau1_phi = tau1_p4.Phi();
  m_fTau1_eta = tau1_p4.Eta();

  m_fTauVis0_pt = tauvis0_p4.Pt() / GeV;
  m_fTauVis0_phi = tauvis0_p4.Phi();
  m_fTauVis0_eta = tauvis0_p4.Eta();
  m_fTauVis1_pt = tauvis1_p4.Pt() / GeV;
  m_fTauVis1_phi = tauvis1_p4.Phi();
  m_fTauVis1_eta = tauvis1_p4.Eta();

  m_fB0_pt = b0_p4.Pt() / GeV;
  m_fB0_phi = b0_p4.Phi();
  m_fB0_eta = b0_p4.Eta();
  m_fB1_pt = b1_p4.Pt() / GeV;
  m_fB1_phi = b1_p4.Phi();
  m_fB1_eta = b1_p4.Eta();

  // deltaRs
  m_fDeltaR_TauTau = tau0_p4.DeltaR(tau1_p4);
  m_fDeltaR_TauVisTauVis = tauvis0_p4.DeltaR(tauvis1_p4);
  m_fDeltaR_BB = b0_p4.DeltaR(b1_p4);
  m_fDeltaR_BB_TauTau = (b0_p4 + b1_p4).DeltaR(tau0_p4 + tau1_p4);

  // Higgs Pt
  m_fPtBB = (b0_p4 + b1_p4).Pt() / GeV;
  m_fPtTauTau = (tau0_p4 + tau1_p4).Pt() / GeV;

  // invariant masses
  m_fMTauTau = (tau0_p4 + tau1_p4).M() / GeV;
  m_fMTauVisTauVis = (tauvis0_p4 + tauvis1_p4).M() / GeV;
  m_fMBB = (b0_p4 + b1_p4).M() / GeV;
  m_fMHH = (tau0_p4 + tau1_p4 + b0_p4 + b1_p4).M() / GeV;

  // if less than two b-jet candidates, just skip this event. If more than 2 b-jet candidates, just take the leading and sub-leading b-jets. You should firstly rank the b-jet in the BjetsJetVec container from largest to smallest according to the pt of the jet
//  if (BjetsJetVec.size() < 2)
//  {
//    return StatusCode::SUCCESS;
//  }
//  else
//  {
//    std::sort(BjetsJetVec.begin(), BjetsJetVec.end(), 
//      [](const xAOD::Jet *a, const xAOD::Jet *b) { return a->pt() > b->pt(); });
//    const xAOD::Jet *bjet0 = BjetsJetVec[0];
//    const xAOD::Jet *bjet1 = BjetsJetVec[1];
//
//    TLorentzVector bjet0_p4, bjet1_p4;
//    bjet0_p4 = bjet0->p4();
//    bjet1_p4 = bjet1->p4();
//    m_fBjet0_pt = bjet0_p4.Pt() / GeV;
//    m_fBjet0_phi = bjet0_p4.Phi();
//    m_fBjet0_eta = bjet0_p4.Eta();
//    m_fBjet1_pt = bjet1_p4.Pt() / GeV;
//    m_fBjet1_phi = bjet1_p4.Phi();
//    m_fBjet1_eta = bjet1_p4.Eta();
//    m_fDeltaR_BjetBjet = bjet0_p4.DeltaR(bjet1_p4);
//    m_fDeltaR_BjetBjet_TauVisTauVis = (bjet0_p4 + bjet1_p4).DeltaR(tauvis0_p4 + tauvis1_p4);
//    m_fMBjetBjet = (bjet0_p4 + bjet1_p4).M() / GeV;
//  }

  // if less than two vbf-tag jet candidates, just skip this event. Otherwise, go to the 2 VBF-tag jets selection block, where 2 methods was used to implemente the selection.
//  if (vbfTagJetVec.size() < 2)
//  {
//    return StatusCode::SUCCESS;
//  }
//  else//take the leading jet and sub-leading jet of the vbfTagJetVec vector as the vbf-tag jets
//  {
//    // firstly, just sort the jets in the vbfTagJetVec from largest to smallest according to the pt of the jet. It's not  necessary at the moment.
//    std::sort(vbfTagJetVec.begin(), vbfTagJetVec.end(), 
//      [](const xAOD::Jet *a, const xAOD::Jet *b) { return a->pt() > b->pt(); });

    // add code for the VBF tag jets with the largest Mjj value (the double for loop and the Mjj sorting out). Idea: to find the jet pair with the largest Mjj.
//l    struct myPoint
//l    { 
//l      double A;
//l      double B;
//l      double C;
//l    };
//l    myPoint mypoint;
//l    std::vector <myPoint> mypoints;
//l    for (std::size_t i = 0; i < vbfTagJetVec.size(); i++)
//l    {
//l      for (std::size_t j = i+1; j < vbfTagJetVec.size(); j++)
//l      {
//l        TLorentzVector tmpVTjet0_p4, tmpVTjet1_p4;
//l        tmpVTjet0_p4 = vbfTagJetVec[i]->p4();
//l        tmpVTjet1_p4 = vbfTagJetVec[j]->p4();
//l        double tmpMVTjetVTjet;
//l        tmpMVTjetVTjet = (tmpVTjet0_p4 + tmpVTjet1_p4).M() / GeV;
        
//l        mypoint.A = i;
//l        mypoint.B = j;
//l        mypoint.C = tmpMVTjetVTjet;
//l        mypoints.push_back(mypoint);
//l      }
//l    }
//l    std::sort(mypoints.begin(), mypoints.end(),
//l      [](const myPoint a, const myPoint b) { return a.C > b.C; });
//l    const xAOD::Jet *VTjet0 = vbfTagJetVec[mypoints[0].A];
//l    const xAOD::Jet *VTjet1 = vbfTagJetVec[mypoints[0].B];

    // add code for the VBF tag jets: 1. group the candidate jets into eta+ and eta-, then choose 2 dedicated jets separately from the 2 groups to make the largest Mjj
//    vector<const xAOD::Jet *> vbfTagJet_etaN_Vec{};
//    vector<const xAOD::Jet *> vbfTagJet_etaP_Vec{};
//    for (std::size_t i = 0; i < vbfTagJetVec.size(); i++)
//    {
//      if ((vbfTagJetVec[i]->p4().Eta()) < 0)
//      {
//        vbfTagJet_etaN_Vec.push_back(vbfTagJetVec[i]);
//      }
//      if ((vbfTagJetVec[i]->p4().Eta()) > 0)
//      { 
//        vbfTagJet_etaP_Vec.push_back(vbfTagJetVec[i]);
//      }
//    }
//    if (vbfTagJet_etaN_Vec.size() < 1 || vbfTagJet_etaP_Vec.size() < 1) return StatusCode::SUCCESS;
//    struct myPoint
//    {
//      double A;
//      double B;
//      double C;
//    };
//    myPoint mypoint;
//    std::vector <myPoint> mypoints;
//    for (std::size_t i = 0; i < vbfTagJet_etaN_Vec.size(); i++)
//    {
//      for (std::size_t j = 0; j < vbfTagJet_etaP_Vec.size(); j++)
//      {
//        TLorentzVector tmpVTjet0_p4, tmpVTjet1_p4;
//        tmpVTjet0_p4 = vbfTagJet_etaN_Vec[i]->p4();
//        tmpVTjet1_p4 = vbfTagJet_etaP_Vec[j]->p4();
//        double tmpMVTjetVTjet;
//        tmpMVTjetVTjet = (tmpVTjet0_p4 + tmpVTjet1_p4).M() / GeV;

//        mypoint.A = i;
//        mypoint.B = j;
//        mypoint.C = tmpMVTjetVTjet;
//        mypoints.push_back(mypoint);
//      }
//    }
//    std::sort(mypoints.begin(), mypoints.end(),
//      [](const myPoint a, const myPoint b) { return a.C > b.C; });
//    const xAOD::Jet *VTjet0 = vbfTagJet_etaN_Vec[mypoints[0].A];
//    const xAOD::Jet *VTjet1 = vbfTagJet_etaP_Vec[mypoints[0].B];
//    ANA_MSG_DEBUG("----------------------------------------is " << VTjet0->p4().Pt());

    //const xAOD::Jet *VTjet0 = vbfTagJetVec[0];
    //const xAOD::Jet *VTjet1 = vbfTagJetVec[1];

//    TLorentzVector VTjet0_p4, VTjet1_p4;
//    VTjet0_p4 = VTjet0->p4();
//    VTjet1_p4 = VTjet1->p4();
//    m_fVTjet0_pt = VTjet0_p4.Pt() / GeV;
//    m_fVTjet0_phi = VTjet0_p4.Phi();
//    m_fVTjet0_eta = VTjet0_p4.Eta();
//    m_fVTjet1_pt = VTjet1_p4.Pt() / GeV;
//    m_fVTjet1_phi = VTjet1_p4.Phi();
//    m_fVTjet1_eta = VTjet1_p4.Eta();
//    m_fDeltaR_VTjetVTjet = VTjet0_p4.DeltaR(VTjet1_p4);
//    m_fDeltaEta_VTjetVTjet = abs(m_fVTjet0_eta - m_fVTjet1_eta);
//    m_fMVTjetVTjet = (VTjet0_p4 + VTjet1_p4).M() / GeV;
//  }

  // select the vbf truth quark
  if (outgoingQuark0Vec.size() < 1 || outgoingQuark1Vec.size() < 1) 
  {
    return StatusCode::SUCCESS;
  }
  else
  {
//    std::sort(outgoingQuark0Vec.begin(), outgoingQuark0Vec.end(),
//      [](const xAOD::TruthParticle *a, const xAOD::TruthParticle *b) { return a->pt() > b->pt(); });
//    std::sort(outgoingQuark1Vec.begin(), outgoingQuark1Vec.end(),
//      [](const xAOD::TruthParticle *c, const xAOD::TruthParticle *d) { return c->pt() > d->pt(); });
//    const xAOD::TruthParticle *VTjet0 = getFinal(outgoingQuark0Vec[0]);
//    const xAOD::TruthParticle *VTjet1 = getFinal(outgoingQuark1Vec[0]);
    const xAOD::TruthParticle *VTjet0 = outgoingQuark0Vec[0];
    const xAOD::TruthParticle *VTjet1 = outgoingQuark1Vec[0];
//    if (VTjet0->p4().Pt() < 10 || VTjet1->p4().Pt() < 10) return StatusCode::SUCCESS;// rule out the events with unreasonable pt value
//    if ((VTjet0->p4().Eta())*(VTjet1->p4().Eta()) >= 0 || abs(VTjet0->p4().Eta()) > 100 || abs(VTjet0->p4().Eta()) < 0.001 || abs(VTjet1->p4().Eta()) > 100 || abs(VTjet1->p4().Eta()) < 0.001) return StatusCode::SUCCESS;// rule out the events with unreasonable eta value
//    ANA_MSG_DEBUG("--------------------------------------Pt is " << VTjet0->p4().Pt());
//    ANA_MSG_DEBUG("--------------------------------------Eta is " << getFinal(outgoingQuark0Vec[0])->p4().Eta());
//    ANA_MSG_DEBUG("--------------------------------------Phi is " << getFinal(outgoingQuark0Vec[0])->p4().Phi());

    TLorentzVector VTjet0_p4, VTjet1_p4;
    VTjet0_p4 = VTjet0->p4();
    VTjet1_p4 = VTjet1->p4();
    m_fVTjet0_pt = VTjet0_p4.Pt() / GeV;
    m_fVTjet0_phi = VTjet0_p4.Phi();
    m_fVTjet0_eta = VTjet0_p4.Eta();
    m_fVTjet1_pt = VTjet1_p4.Pt() / GeV;
    m_fVTjet1_phi = VTjet1_p4.Phi();
    m_fVTjet1_eta = VTjet1_p4.Eta();
    m_fDeltaR_VTjetVTjet = VTjet0_p4.DeltaR(VTjet1_p4);
    m_fDeltaEta_VTjetVTjet = abs(m_fVTjet0_eta - m_fVTjet1_eta);
    m_fMVTjetVTjet = (VTjet0_p4 + VTjet1_p4).M() / GeV;
  }

  // apply addition cuts for VBF tag jets
//l  APPLYCUT(m_fMVTjetVTjet > 1000, "Mjj cut");
//l  APPLYCUT(m_fDeltaEta_VTjetVTjet > 3, "DeltaEta_jj");
//l  APPLYCUT(m_fVTjet0_eta * m_fVTjet1_eta < 0, "Eta*EtaOS_jj");

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
