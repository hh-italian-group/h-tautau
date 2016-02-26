/*!
 * \file ElectronsIDAnalyzer.cc
 * \author Original author: Ilya Kravchenko
 * \author Contributing author: Claudio Caputo
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HHbbTauTau/AnalysisBase/include/SyncTree.h"


#include "TTree.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class ElectronsIDAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ElectronsIDAnalyzer(const edm::ParameterSet&);
      ~ElectronsIDAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0,
              TRUE_PROMPT_ELECTRON,
              TRUE_ELECTRON_FROM_TAU,
              TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,
               const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
                    int &ancestorPID, int &ancestorStatus);

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      // ... none ...

      bool electronBool_;
      // AOD case data members
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

      //Tau Tag
      edm::EDGetToken tausMiniAODToken_;

  TTree *electronTree_;
  TTree *tauTree_;

  // Global info
  Int_t run_;
  Int_t lumi_;
  Int_t evtnum_;

  // all variables for the output tree

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  //Electrons
  Int_t nElectrons_;

  std::vector<Float_t> mvaValue_;
  std::vector<Int_t>   mvaCategory_;

  std::vector<Int_t> passMediumId_;
  std::vector<Int_t> passTightId_;

  std::vector<Int_t> isTrue_;

  //Taus
  Int_t nTau_;

  std::vector<Int_t> decayMode_;
  std::vector<Int_t> decayModeFindingOld_;
  std::vector<Int_t> decayModeFindingNew_;

  std::vector<Float_t> isoCombinedRaw_;
  std::vector<Float_t> isoMVARawOldDM_;
  std::vector<Float_t> isoMVARawNewDM_;

  std::vector<Int_t> muonLoose_;
  std::vector<Int_t> muonTight_;

  std::vector<Int_t> electronVLoose_;
  std::vector<Int_t> electronLoose_;
  std::vector<Int_t> electronMedium_;

  ntuple::SyncTree syncTree;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronsIDAnalyzer::ElectronsIDAnalyzer(const edm::ParameterSet& iConfig):
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))),
  syncTree(&edm::Service<TFileService>()->file(),false)
{

  //
  // Prepare tokens for all input collections and objects
  //
  electronBool_ = iConfig.getUntrackedParameter<bool>("electronBool");
  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into
  // GsfElectron
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

    tausMiniAODToken_= mayConsume<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"));


  //edm::Service<TFileService> fs;

  if(electronBool_){
//  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");


//  electronTree_->Branch("run"        ,  &run_     , "run/I");
//  electronTree_->Branch("lumi"       ,  &lumi_    , "lumi/I");
//  electronTree_->Branch("evtnum"     ,  &evtnum_  , "evtnum/I");

//  electronTree_->Branch("nEle",  &nElectrons_ , "nEle/I");
//  electronTree_->Branch("pt"  ,  &pt_    );
//  electronTree_->Branch("eta" ,  &eta_ );
//  electronTree_->Branch("phi" ,  &phi_ );

//  electronTree_->Branch("mvaVal" ,  &mvaValue_ );
//  electronTree_->Branch("mvaCat" ,  &mvaCategory_ );

//  electronTree_->Branch("passMediumId" ,  &passMediumId_ );
//  electronTree_->Branch("passTightId"  ,  &passTightId_ );

//  electronTree_->Branch("isTrue"             , &isTrue_);
  }
  // Tau TTree
  if(!electronBool_){
//  tauTree_      = fs->make<TTree> ("TauTree", "Tau data");

//  tauTree_->Branch("run"        ,  &run_     , "run/I");
//  tauTree_->Branch("lumi"       ,  &lumi_    , "lumi/I");
//  tauTree_->Branch("evtnum"     ,  &evtnum_  , "evtnum/I");

//  tauTree_->Branch("nTau",  &nTau_ , "nTau/I");
////  tauTree_->Branch("pt"  ,  &pt_    );
////  tauTree_->Branch("eta" ,  &eta_ );
////  tauTree_->Branch("phi" ,  &phi_ );

//  tauTree_->Branch("decayMode" ,  &decayMode_ );
//  tauTree_->Branch("decayModeFindingOldDMs" ,  &decayModeFindingOld_ );
//  tauTree_->Branch("decayModeFindingNewDMs" ,  &decayModeFindingNew_ );

//  tauTree_->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits" ,  &isoCombinedRaw_ );
//  tauTree_->Branch("byIsolationMVA3oldDMwLTraw" ,  &isoMVARawOldDM_ );
//  tauTree_->Branch("byIsolationMVA3newDMwLTraw" ,  &isoMVARawNewDM_ );

//  tauTree_->Branch("againstMuonLoose3"  ,  &muonLoose_ );
//  tauTree_->Branch("againstMuonTight3"  ,  &muonTight_ );

//  tauTree_->Branch("againstElectronVLooseMVA5"  ,  &electronVLoose_ );
//  tauTree_->Branch("againstElectronLooseMVA5"  ,  &electronLoose_ );
//  tauTree_->Branch("againstElectronMediumMVA5"  ,  &electronMedium_ );
  }
}


ElectronsIDAnalyzer::~ElectronsIDAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronsIDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Save global info right away
  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  evtnum_ = iEvent.id().event();

  syncTree.run() = iEvent.id().run();
  syncTree.evt() = iEvent.id().event();

  // Retrieve the collection of electrons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name.
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }

  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  //Get Tau collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tausMiniAODToken_, taus);
  //const pat::TauCollection& taus = *hTaus.product();

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  if (!electronBool_){// Clear Taus vectors
  nTau_ = 0;
//  pt_.clear();
//  eta_.clear();
//  phi_.clear();
  //
  decayMode_.clear();
  decayModeFindingOld_.clear();
  decayModeFindingNew_.clear();

  isoCombinedRaw_.clear();
  isoMVARawNewDM_.clear();
  isoMVARawOldDM_.clear();

  muonLoose_.clear();
  muonTight_.clear();
  electronVLoose_.clear();
  electronLoose_.clear();
  electronMedium_.clear();

  //loop over taus
  //for (const auto tau : taus){
  for(size_t i = 0; i < taus->size(); ++i){
      const auto tau = taus->ptrAt(i);
      if(tau->pt() > 20 && std::abs(tau->eta()) < 2.3 ) {
          nTau_++;

//          pt_   .push_back(tau->pt());
//          eta_  .push_back(tau->eta());
//          phi_  .push_back(tau->phi());
          decayMode_ .push_back(tau->decayMode());

          decayModeFindingOld_.push_back(tau->tauID("decayModeFinding"));
          decayModeFindingNew_.push_back(tau->tauID("decayModeFindingNewDMs"));

          isoCombinedRaw_     .push_back(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
          isoMVARawOldDM_     .push_back(tau->tauID("byIsolationMVA3oldDMwLTraw"));
          isoMVARawOldDM_     .push_back(tau->tauID("byIsolationMVA3newDMwLTraw"));


          muonLoose_          .push_back(tau->tauID("againstMuonLoose3"));
          muonTight_          .push_back(tau->tauID("againstMuonTight3"));
          electronVLoose_     .push_back(tau->tauID("againstElectronVLooseMVA5"));
          electronLoose_      .push_back(tau->tauID("againstElectronLooseMVA5"));
          electronMedium_     .push_back(tau->tauID("againstElectronMediumMVA5"));


      }
    if(nTau_ == 1){
        syncTree.pt_1() = tau->pt();
        syncTree.eta_1() = tau->eta();
        syncTree.phi_1() = tau->phi();
    }
  }
  //tauTree_->Fill();
  }

  if(electronBool_){
  // Clear Electron vectors
  nElectrons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  //
  mvaValue_.clear();
  mvaCategory_.clear();
  passMediumId_.clear();
  passTightId_ .clear();
  //
  isTrue_.clear();

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics
    if( el->pt() > 20 && std::abs(el->eta()) < 2.0 ) { // keep only electrons above 5 GeV

     nElectrons_++;

    //
    // Save electron kinematics
    //
     pt_  .push_back( el->pt() );
     eta_ .push_back( el->superCluster()->eta() );
     phi_ .push_back( el->superCluster()->phi() );

    //
    // Look up and save the ID decisions
    //
     bool isPassMedium = (*medium_id_decisions)[el];
     bool isPassTight  = (*tight_id_decisions)[el];
     passMediumId_.push_back( (int)isPassMedium);
     passTightId_.push_back ( (int)isPassTight );

     mvaValue_.push_back( (*mvaValues)[el] );
     mvaCategory_.push_back( (*mvaCategories)[el] );

    // Save MC truth match
     isTrue_.push_back( matchToTruth( el, genParticles) );

    }
   }

  // Save the info
  //electronTree_->Fill();
  }

  syncTree.Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
ElectronsIDAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronsIDAnalyzer::endJob()
{
    syncTree.Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
ElectronsIDAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronsIDAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronsIDAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronsIDAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronsIDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int ElectronsIDAnalyzer::matchToTruth(const edm::Ptr<reco::GsfElectron> el,
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  //
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  //
  int ancestorPID = -999;
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }

  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronsIDAnalyzer::findFirstNonElectronMother(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronsIDAnalyzer);
