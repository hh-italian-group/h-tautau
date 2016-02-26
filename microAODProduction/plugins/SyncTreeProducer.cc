/*!
 * \file SyncTreeProducer.cc
 * \author author: Claudio Caputo
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


// import LHEEventProduction definition
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//SVFit
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

//HHbbTauTau Framework
#include "../../TreeProduction/interface/Tau.h"
#include "../../TreeProduction/interface/Muon.h"
#include "../../AnalysisBase/include/AnalyzerData.h"
#include "../../AnalysisBase/include/CutTools.h"
#include "../interface/Candidate.h"
#include "../interface/SyncTree.h"
#include "../../Analysis/include/SelectionResults.h"

#include "../interface/Htautau_2015.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

using namespace analysis;

//Analyzer Data Class
namespace analysis {
class SyncAnalyzerData : public root_ext::AnalyzerData {
public:
    SyncAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)

};

struct SelectionResultsV2_mutau : public SelectionResultsV2 {
    finalState::bbMuTaujet muTau_MC;
    CandidateV2Ptr GetMuon() const { return higgs->GetDaughter(CandidateV2::Type::Muon); }
    CandidateV2Ptr GetTau() const { return higgs->GetDaughter(CandidateV2::Type::Tau); }

    virtual CandidateV2Ptr GetLeg(size_t leg_id) const override
    {
        if(leg_id == 1) return GetMuon();
        if(leg_id == 2) return GetTau();
        throw exception("Bad leg id = ") << leg_id;
    }

    //virtual const finalState::bbTauTau& GetFinalStateMC() const override { return muTau_MC; }

};

}
//
// class declaration
//

class SyncTreeProducer : public edm::EDAnalyzer {
   public:
      explicit SyncTreeProducer(const edm::ParameterSet&);
      ~SyncTreeProducer();

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

      analysis::SyncAnalyzerData& GetAnaData() { return anaData; }

      void FillSyncTree(const edm::Event&);

      analysis::CandidateV2PtrVector FindCompatibleObjects(const CandidateV2PtrVector& objects1,
                                                           const CandidateV2PtrVector& objects2,
                                                           double minDeltaR, CandidateV2::Type type,
                                                           const std::string& hist_name,
                                                           int expectedCharge=analysis::CandidateV2::UnknownCharge());

      analysis::CandidateV2PtrVector ApplyTriggerMatch(pat::TriggerObjectStandAloneCollection triggerObjects,
                                                       const CandidateV2PtrVector& higgses,const edm::TriggerNames &names,
                                                       const std::set<std::string>& hltPaths,
                                                       bool useStandardTriggerMatch, const bool isCrossTrigger);


      double Fit(const CandidateV2Ptr& higgs, const analysis::MissingETPtr& met);
      analysis::CandidateV2Ptr SelectSemiLeptonicHiggs(CandidateV2PtrVector& higgses);

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,
               const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
                    int &ancestorPID, int &ancestorStatus);



      inline bool HaveTriggerMatched(pat::TriggerObjectStandAloneCollection& triggerObjects,
                                     const edm::TriggerNames &names,
                                     const std::set<std::string>& interestingPaths, const CandidateV2& candidate,
                                     const double deltaR_Limit,const bool isCrossTrigger)
      {
          for (const std::string& interestinPath : interestingPaths){
              std::cout<<"Trigger Path  :   "<<interestinPath<<std::endl;
              if (HaveTriggerMatched(triggerObjects,names,interestinPath,candidate, deltaR_Limit, isCrossTrigger)) return true;
          }
          return false;
      }
//Now it could manage also the single object Trigger
      inline bool HaveTriggerMatched(pat::TriggerObjectStandAloneCollection& triggerObjects,
                                     const edm::TriggerNames &names,
                                     const std::string& interestingPath,
                                     const CandidateV2& candidate, const double deltaR_Limit, const bool isCrossTrigger)
      {
          if(candidate.GetFinalStateDaughters().size()) {
              for(const auto& daughter : candidate.GetFinalStateDaughters()) {
                  if(isCrossTrigger &&
                          !SyncTreeProducer::HaveTriggerMatched(triggerObjects, names, interestingPath, *daughter, deltaR_Limit,isCrossTrigger))
                      return false;
                  if(!isCrossTrigger &&
                          daughter->GetType() == CandidateV2::Type::Muon &&
                          SyncTreeProducer::HaveTriggerMatched(triggerObjects, names, interestingPath, *daughter, deltaR_Limit,isCrossTrigger))
                      return true;
              }
              if ( isCrossTrigger) return true;
              if (!isCrossTrigger) return false;
          }

          for (auto &triggerObject : triggerObjects){
              triggerObject.unpackPathNames(names);
              TLorentzVector triggerObjectMomentum;
              triggerObjectMomentum.SetPtEtaPhiM(triggerObject.pt(), triggerObject.eta(), triggerObject.phi(), triggerObject.mass());
              //Pt cut for singleTrigger
              // if (!isCrossTrigger && triggerObject.pt() < 18) continue;

              for (unsigned n = 0; n < triggerObject.pathNames(true).size(); ++n){
                  const std::string& objectMatchedPath = triggerObject.pathNames(true).at(n);

                  const size_t found = objectMatchedPath.find(interestingPath);
                  bool isBoth = triggerObject.hasPathName( triggerObject.pathNames(true).at(n), true, true );
                  if (found != std::string::npos) {
                      std::cout<<"\t Trigger obj :   "<<triggerObjectMomentum<<std::endl;
                      std::cout << "\t\t Path Names :   " << objectMatchedPath
                            << (isBoth ? "  PASS\n " : "  fail (or not run) --\n ")
                            //<< triggerObject.pathLastFilterAccepted.at(n) << "\n"
                            << "\t\t Candidato :  " << candidate.GetMomentum()
                            <<"\t\t DeltaR = "<< triggerObjectMomentum.DeltaR(candidate.GetMomentum())
                            <<(triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit ? "\t Matched":"\t NotMatched")
                            << std::endl;
                  }
                  if (found != std::string::npos && isBoth &&
                          triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit)
                      return true;
              }
          }
          return false;
      }

    //  bool HaveTriggerFired(const std::set<std::string>& interestinghltPaths) const ;

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      // ... none ...

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
  //    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

      //InputTag
      edm::EDGetToken tausMiniAODToken_;
      edm::EDGetToken muonsMiniAODToken_;
      edm::EDGetToken vtxMiniAODToken_;
      edm::EDGetToken pfMETAODToken_;
      edm::EDGetToken jetsMiniAODToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUInfo_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrixTAG_;
      edm::EDGetTokenT<LHEEventProduct> lheEventProduct_;
      edm::EDGetTokenT<GenEventInfoProduct> genWeights_;

      bool computeHT_;
      Run2::SyncTree syncTree;
      analysis::SyncAnalyzerData anaData;
      SelectionResultsV2_mutau selection;
      VertexV2Ptr primaryVertex;
      std::string sampleType;
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
SyncTreeProducer::SyncTreeProducer(const edm::ParameterSet& iConfig):
  tausMiniAODToken_(mayConsume<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
  muonsMiniAODToken_(mayConsume<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  vtxMiniAODToken_(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
  pfMETAODToken_(mayConsume<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
  jetsMiniAODToken_(mayConsume<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  PUInfo_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  metCovMatrixTAG_(consumes
                   <ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >>(iConfig.getParameter<edm::InputTag>("metCov"))),
  lheEventProduct_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
  genWeights_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  computeHT_(iConfig.getParameter<bool>("HTBinning")),
  syncTree(&edm::Service<TFileService>()->file(),false),
  anaData("BeforCut.root"),
  sampleType(iConfig.getParameter<std::string>("sampleType"))
{
//  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
//    (iConfig.getParameter<edm::InputTag>
//     ("genParticlesMiniAOD"));

}


SyncTreeProducer::~SyncTreeProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::MuTau;

  cuts::Cutter cut(&GetAnaData().Selection("events"));

  // Save global info right away
  //Get collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tausMiniAODToken_, taus);
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsMiniAODToken_, muons);
  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(vtxMiniAODToken_, vertices);
  edm::Handle<edm::View<pat::MET> > pfMETs;
  iEvent.getByToken(pfMETAODToken_, pfMETs);
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetsMiniAODToken_, jets);
  edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrix;
  iEvent.getByToken(metCovMatrixTAG_,metCovMatrix);
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
  iEvent.getByToken(PUInfo_, PUInfo);

  try{

      std::cout<< "=========================================================================== \n"
               << "\t Run "<<iEvent.id().run()<<"\t Lumi "<<iEvent.id().luminosityBlock()
                  << "\t Event "<<iEvent.id().event()<<std::endl;

      cut(true,"events");

      if(computeHT_){
          //access generator level HT
      edm::Handle<LHEEventProduct> lheEventProduct;
      iEvent.getByToken(lheEventProduct_, lheEventProduct);
      const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
      double lheHt = 0.;
      size_t numParticles = lheParticles.size();
      for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
       int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
       int status = lheEvent.ISTUP[idxParticle];
       if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
           lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
       }
      }
      selection.HT = lheHt;
      //std::cout<<"  HT  ============ > "<<lheHt<<std::endl;
      }

      if (sampleType == "Spring15MC" || sampleType == "Fall15MC"){
          edm::Handle<GenEventInfoProduct> genEvt;
          iEvent.getByToken(genWeights_,genEvt);
          // event weight
          selection.weightevt=genEvt->weight();
      }


          if(PUInfo.isValid()){
              for(std::vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI){
                  int BX = PVI->getBunchCrossing();
//                  std::cout << "\t PileUpInfo :    BX = " << PVI->getBunchCrossing()
//                            << "  PU = " << PVI->getTrueNumInteractions() << std::endl;
                  if(BX == 0) selection.numtruepileupinteractions = PVI->getTrueNumInteractions();
            }
          }

      bool triggerFired = false;
      const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
      const auto Key = analysis::stringToDataSourceTypeMap.at(sampleType);
      const auto& hltPaths = MuTau::trigger::hltPathMaps.at(Key);
    //  const auto& hltPaths = MuTau::trigger::hltPathMC;
          //std::cout << "\n === TRIGGER PATHS === " << std::endl;
          for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
              for (const std::string& triggerPath : hltPaths ){

                  const std::string& objectMatchedPath = names.triggerName(i);
                  size_t found = objectMatchedPath.find(triggerPath);


                  if(found != std::string::npos){
//                      std::cout << "TriggerPath --->   " << triggerPath << std::endl;
//                      std::cout << "Trigger " << names.triggerName(i) <<
//                      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
//                      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
//                      << std::endl;
                      if(triggerPrescales->getPrescaleForIndex(i) == 1 && triggerBits->accept(i)) triggerFired = true;
                }
              }
          }

      cut(triggerFired,"trigger");
  
      const auto PV = (*vertices).ptrAt(0); //Deferenzio per passare da edm::Handle a edm::View. Quest'ultimo permette
                                            //di gestire una qualsiasi collezione del tipo passatogli tramite Tamplate.
                                            //Es. edm::View<int> gestisce int semplici, vector<int>, set<int> etc.

     // auto primVtx = *(PV.get());
      for (const reco::Vertex& vertex: *vertices){

//          if (!(vertex.ndof() < cuts::Htautau_2015::vertex::ndf &&
//                std::abs(vertex.z()) < cuts::Htautau_2015::vertex::z )) continue;
          VertexV2Ptr tmp_Vertex(new VertexV2(vertex));
          selection.vertices.push_back(tmp_Vertex);
      }

      cut(selection.vertices.size(),"vertex");

      primaryVertex = selection.GetPrimaryVertex();
      std::cout<<" Vertici    --->   "<<selection.vertices.size()<<std::endl;

      std::vector<pat::Muon> patMuonsVector;
      CandidateV2PtrVector muonCollection, tauCollection;
      CandidateV2PtrVector muonVetoCollection;
      for(const pat::Muon &muon : *muons){
        if(!(muon.pt() > muonID::pt  &&
            fabs(muon.eta()) < muonID::eta &&
            muon.isMediumMuon() == muonID::isMediumMuon)) continue;

        bool muonIP = fabs(muon.muonBestTrack()->dxy(PV->position())) < muonID::dB &&
                      fabs(muon.muonBestTrack()->dz(PV->position())) < muonID::dz;
        if(!muonIP) continue;

        const CandidateV2Ptr muon_candidate(new CandidateV2(muon));
        muonCollection.push_back(muon_candidate);
        patMuonsVector.push_back(muon);

      }

      for(const pat::Muon &muon : *muons){
          double iso_mu = (muon.pfIsolationR03().sumChargedHadronPt + std::max(
                            muon.pfIsolationR03().sumNeutralHadronEt +
                            muon.pfIsolationR03().sumPhotonEt -
                            0.5 * muon.pfIsolationR03().sumPUPt, 0.0)) / muon.pt();
          bool muonIP = fabs(muon.muonBestTrack()->dxy(PV->position())) < muonID::dB &&
                        fabs(muon.muonBestTrack()->dz(PV->position())) < muonID::dz;
          if(!(muon.pt() > ZmumuVeto::pt && fabs(muon.eta()) < ZmumuVeto::eta &&
               muon.isGlobalMuon() && muon.isTrackerMuon()  && muon.isPFMuon()  &&
               muonIP && iso_mu < 0.3)) continue;

          const CandidateV2Ptr muon_candidate(new CandidateV2(muon));
          muonVetoCollection.push_back(muon_candidate);
      }

      cut(muonCollection.size(),"muons");

      auto Zmumu = FindCompatibleObjects(muonVetoCollection,muonVetoCollection,DeltaR_betweenSignalObjects,
                                         CandidateV2::Type::Z,"Z_mu_mu",0);
      selection.Zveto = Zmumu.size() ? true : false;

      std::cout<< "Muons Checks ------------------------------------------------ \n";
      for(auto &muon : muonCollection){
          const pat::Muon& patMuon = muon->GetNtupleObject<pat::Muon>();
          double iso_mu = (patMuon.pfIsolationR03().sumChargedHadronPt + std::max(
                            patMuon.pfIsolationR03().sumNeutralHadronEt +
                            patMuon.pfIsolationR03().sumPhotonEt -
                            0.5 * patMuon.pfIsolationR03().sumPUPt, 0.0)) / patMuon.pt();
          std::cout << "\t Candidate Muon :  "<< muon->GetMomentum()<<"\n"
                       << "\t Pat Muon  Casted  "<<&patMuon<<" :  "<< patMuon.pt() << "  "<< patMuon.eta()
                       << "  Iso : "<< iso_mu<<"\n";

      }

      //Taus selections
      for (const pat::Tau &tau : *taus){

      pat::PackedCandidate const* packedLeadTauCand =
              dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      //fabs(packedLeadTauCand->dz()) < 0.2;  // The PackedCandidate::dz() method is wrt. the first PV by default


          if(!(tau.pt() > tauID::pt && fabs(tau.eta()) < tauID::eta &&
               tau.tauID("decayModeFinding") > tauID::decayModeFinding &&
               fabs(packedLeadTauCand->dz()) < tauID::dz )) continue;

	  if(!(fabs(tau.charge())==1)) continue;

        const CandidateV2Ptr tau_candidate(new CandidateV2(tau));

        tauCollection.push_back(tau_candidate);
      }

    cut(tauCollection.size(),"taus");

    std::cout<< "Taus Checks ------------------------------------------------ \n";
    for(auto &tau : tauCollection){
        const pat::Tau& patTau = tau->GetNtupleObject<pat::Tau>();
        double iso_tau = patTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        std::cout << "\t Candidate Tau :  "<< tau->GetMomentum()<<"\n"
                     << "\t Pat Tau  Casted  "<<&patTau<<" : "<< patTau.pt() << "  "<< patTau.eta()
                     << "  Iso : "<< iso_tau<<"\n";

    }


    auto higgses = FindCompatibleObjects(muonCollection,tauCollection,DeltaR_betweenSignalObjects,CandidateV2::Type::Higgs,
                                         "H_mu_tau");

    cut(higgses.size(),"mu_tau");

    bool isCrossTrigger = false;
    if ( sampleType == "Run2015C" ) isCrossTrigger = true;

    auto triggeredHiggses = ApplyTriggerMatch(*(triggerObjects.product()), higgses,names,hltPaths,false,isCrossTrigger);

    cut(triggeredHiggses.size(),"triggerMatch");

    for (auto hig : triggeredHiggses){
        std::cout<<" Triggered Pair : "<< hig->GetMomentum() <<std::endl;
    }

    selection.higgs = SelectSemiLeptonicHiggs(triggeredHiggses);
    auto higgs = selection.higgs;


    std::cout<< "\n\t CHOOSEN -->  " << higgs->GetMomentum() <<std::endl;
//    const pat::Tau&  tauBuono  = higgs->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
//    std::cout<< "\n\t Tau Iso-->  " << tauBuono.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") <<std::endl;


    MissingETPtr pfMET(new MissingET((*pfMETs)[0],*metCovMatrix));
    selection.pfMET = pfMET;


    CandidateV2PtrVector jetsCollection, bjetsCollection;
    for(const pat::Jet &jet : *jets){

        if (!( jet.pt()>jetID::pt_loose && fabs(jet.eta())<jetID::eta && jetID::passPFLooseId(jet) ) ) continue;

        const CandidateV2Ptr jet_candidate(new CandidateV2(jet));

        if (!(jet_candidate->GetMomentum().DeltaR(selection.GetLeg(1)->GetMomentum()) > jetID::deltaR_signalObjects &&
            jet_candidate->GetMomentum().DeltaR(selection.GetLeg(2)->GetMomentum()) > jetID::deltaR_signalObjects) ) continue;
        jetsCollection.push_back(jet_candidate);

        if ( jet.pt()>jetID::pt_loose && fabs(jet.eta())<btag::eta &&
             jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btag::CSV)
            bjetsCollection.push_back(jet_candidate);
    }

   selection.bjets = bjetsCollection;
   selection.jets  = jetsCollection;

    double svfit_mass = Fit(higgs,pfMET);
    std::cout<< "\n\t CHOOSEN SVFit -->  " << svfit_mass <<std::endl;
    std::cout<< "=========================================================================== \n";


    FillSyncTree(iEvent);
  }catch(cuts::cut_failed&){}
  GetAnaData().Selection("events").fill_selection();
  selection.vertices.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void
SyncTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer::endJob()
{
    syncTree.Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
SyncTreeProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
SyncTreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SyncTreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SyncTreeProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SyncTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SyncTreeProducer::matchToTruth(const edm::Ptr<reco::GsfElectron> el,
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

void SyncTreeProducer::findFirstNonElectronMother(const reco::Candidate *particle,
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

///////////


analysis::CandidateV2PtrVector SyncTreeProducer::FindCompatibleObjects(const CandidateV2PtrVector& objects1,
                                                   const CandidateV2PtrVector& objects2,
                                                   double minDeltaR, CandidateV2::Type type, const std::string& hist_name,
                                                   int expectedCharge)
    {
        CandidateV2PtrVector result;
        for(const CandidateV2Ptr& object1 : objects1) {
            for(const CandidateV2Ptr& object2 : objects2) {
                std::cout << " \t\t\t  Pair DeltaR = "
                          << object2->GetMomentum().DeltaR(object1->GetMomentum()) << std::endl;
                if(object2->GetMomentum().DeltaR(object1->GetMomentum()) > minDeltaR) {
                    const CandidateV2Ptr candidate(new CandidateV2(type, object1, object2));
                    if (expectedCharge != CandidateV2::UnknownCharge() && candidate->GetCharge() != expectedCharge)
                        continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),1);
                }
            }
        }

        GetAnaData().N_objects(hist_name).Fill(result.size(),1);
        return result;
    }

CandidateV2PtrVector SyncTreeProducer::ApplyTriggerMatch(pat::TriggerObjectStandAloneCollection triggerObjects,
                                                         const CandidateV2PtrVector& higgses, const edm::TriggerNames &names,
                                                         const std::set<std::string>& hltPaths,
                                                         bool useStandardTriggerMatch, const bool isCrossTrigger)
{
    CandidateV2PtrVector triggeredHiggses;
    for (const auto& higgs : higgses){

        std::cout<<"### Higgs Pair :  \n"<<higgs->GetMomentum()<<std::endl;
        if(!useStandardTriggerMatch && SyncTreeProducer::HaveTriggerMatched(triggerObjects,names, hltPaths, *higgs,
                                                                    cuts::Htautau_2015::DeltaR_triggerMatch, isCrossTrigger)){
            std::cout<<"### Pushed Back "<<std::endl;
            triggeredHiggses.push_back(higgs);
        }
//        if (useStandardTriggerMatch && SyncTreeProducer::HaveTriggerMatched(hltPaths, *higgs))
//            triggeredHiggses.push_back(higgs);
    }
    return triggeredHiggses;
}

CandidateV2Ptr SyncTreeProducer::SelectSemiLeptonicHiggs(CandidateV2PtrVector& higgses)
{
    if(!higgses.size())
        throw std::runtime_error("no available higgs candidate to select");
    if(higgses.size()==1) return higgses.front();

    const auto higgsSelector = [&] (const CandidateV2Ptr& first, const CandidateV2Ptr& second) -> bool
    {
        const pat::Muon& first_muon = first->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  first_tau  = first->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
        const pat::Muon& second_muon = second->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  second_tau  = second->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();

        double iso_mu1 = (first_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          first_muon.pfIsolationR03().sumNeutralHadronEt +
                          first_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * first_muon.pfIsolationR03().sumPUPt, 0.0)) / first_muon.pt();
        double iso_mu2 = (second_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          second_muon.pfIsolationR03().sumNeutralHadronEt +
                          second_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * second_muon.pfIsolationR03().sumPUPt, 0.0)) / second_muon.pt();

        double iso_tau1 = first_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        double iso_tau2 = second_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");

        std::cout << "LAMBDA ------------------------  \n IsoMu1 = "<<iso_mu1<<"  IsoMu2 = "<<iso_mu2
                  << "  IsoTau1 = "<<iso_tau1<<"  IsoTau2 = "<<iso_tau2<<std::endl;
        std::cout << " PtMu1 = "<<first_muon.pt()<<"  PtMu2 = "<<second_muon.pt()
                  << "  PtTau1 = "<<first_tau.pt()<<"  PtTau2 = "<<second_tau.pt()<<std::endl;

        if (!iso_mu1 || !iso_mu2 ){
            std::cout << "\t 1st Muon \t ChargedHadronPt = "<<first_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<first_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<first_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<first_muon.pfIsolationR03().sumPUPt<<std::endl;
            std::cout << "\t 2nd Muon \t ChargedHadronPt = "<<second_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<second_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<second_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<second_muon.pfIsolationR03().sumPUPt<<std::endl;
        }
        bool muon = (iso_mu1 < iso_mu2) ||
                    ((iso_mu1 == iso_mu2) ? first_muon.pt() >= second_muon.pt() : false);
        bool tau = (iso_tau1 > iso_tau2) ||
                    ((iso_tau1 == iso_tau2) ? first_tau.pt() >= second_tau.pt() : false);

//        if ( &first_muon==&second_muon ){
//            if ( &first_tau==&second_tau ) return true;
//            return tau;
//        }

//        if (!muon) return tau;
//        return muon;
        if ( &first_muon!=&second_muon ) return muon;
        if ( &first_tau!=&second_tau ) return tau;

        throw analysis::exception("not found a good criteria for best tau pair");
    };

    std::sort(higgses.begin(), higgses.end(), higgsSelector) ;
    return higgses.front();
}

double SyncTreeProducer::Fit(const CandidateV2Ptr& higgs, const analysis::MissingETPtr& met){

    const pat::Muon& muon = higgs->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
    const pat::Tau&  tau  = higgs->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
    // define MET
      double measuredMETx = met->Px();
      double measuredMETy = met->Py();
      // define MET covariance
      TMatrixD covMET(2, 2);
      covMET[0][0] = met->GetCovVector().at(0);
      covMET[1][0] = met->GetCovVector().at(1);
      covMET[0][1] = met->GetCovVector().at(2);
      covMET[1][1] = met->GetCovVector().at(3);
      // define lepton four vectors
      std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, muon.pt(),
                                                                      muon.eta(), muon.phi(), muon.mass())); // tau -> electron decay (Pt, eta, phi, mass)
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, tau.pt(),
                                                                      tau.eta(), tau.phi(), tau.mass(),tau.decayMode())); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
      SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
      algo.addLogM(false);
      edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
      TH1::AddDirectory(false);
      TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
      algo.shiftVisPt(true, inputFile_visPtResolution);
      algo.integrateMarkovChain();

      double mass = algo.getMass(); // return value is in units of GeV
      if ( algo.isValidSolution() ) {
          selection.svfitResult.mass = algo.mass();
          selection.svfitResult.has_valid_mass = true;
          //if(fitAlgorithm == FitAlgorithm::MarkovChain) {
              selection.svfitResult.momentum.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.mass());
              selection.svfitResult.has_valid_momentum = true;
         // }
         std::cout << "... m svfit : " << algo.mass() << " +/- " << algo.massUncert() << std::endl;
      } else {
        std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
      }

      delete inputFile_visPtResolution;

      return mass;
}

void SyncTreeProducer::FillSyncTree(const edm::Event& iEvent)
    {
        static const float default_value = Run2::DefaultFloatFillValueForSyncTree();

        // Event
        std::cout<<"~~~~~~~~~~~~~EVENT Info~~~~~~~~~"<<std::endl;
        std::cout<<"Run = "<<iEvent.id().run()<<"  Lumi = "<<iEvent.id().luminosityBlock()
                <<" Event = "<<iEvent.id().event()<<std::endl;
        syncTree.run()  = iEvent.id().run();
        syncTree.lumi() = iEvent.id().luminosityBlock();
        syncTree.evt()  = iEvent.id().event();

        if(computeHT_){
            syncTree.HT()   = selection.HT;
            if(selection.HT<100) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::lt100);
            if(100<=selection.HT && selection.HT<200) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f100to200);
            if(200<=selection.HT && selection.HT<400) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f200to400);
            if(400<=selection.HT && selection.HT<600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f400to600);
            if(selection.HT>=600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::gt600);
//        syncTree->eventType() = static_cast<int>(selection.eventType);
//        syncTree->eventEnergyScale() = static_cast<int>(eventEnergyScale);
        }
        else {
            syncTree.HT() = default_value;
            syncTree.HTBin() = -1;
        }

        if (sampleType == "Spring15MC" || sampleType == "Fall15MC") syncTree.weightevt() = selection.weightevt;
        else syncTree.weightevt() = default_value;

        syncTree.npv() = selection.vertices.size();
        syncTree.npu() = selection.numtruepileupinteractions;
//        if (config.ApplyPUreweight()){
//            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
//            if(bxIndex >= event->eventInfo().bunchCrossing.size())
//                throw std::runtime_error("in-time BX not found");
//            syncTree->npu() = event->eventInfo().trueNInt.at(bxIndex);
//        }

        // Weights
//        syncTree->puweight() = GetEventWeights().GetPileUpWeight();
//        syncTree->trigweight_1() = GetEventWeights().GetTriggerWeight(1);
//        syncTree->trigweight_2() = GetEventWeights().GetTriggerWeight(2);
//        syncTree->idweight_1() = GetEventWeights().GetIdWeight(1);
//        syncTree->idweight_2() = GetEventWeights().GetIdWeight(2);
//        syncTree->isoweight_1() = GetEventWeights().GetIsoWeight(1);
//        syncTree->isoweight_2() = GetEventWeights().GetIsoWeight(2);
//        syncTree->fakeweight_1() = GetEventWeights().GetFakeWeight(1); //
//        syncTree->fakeweight_2() = GetEventWeights().GetFakeWeight(2); // e -> tau fake rate * jet -> tau fake rate
//        syncTree->weight() = GetEventWeights().GetFullWeight();
//        syncTree->embeddedWeight() = GetEventWeights().GetEmbeddedWeight();
//        syncTree->decayModeWeight_1() = GetEventWeights().GetDecayModeWeight(1);
//        syncTree->decayModeWeight_2() = GetEventWeights().GetDecayModeWeight(2);

        // HTT candidate
        syncTree.m_vis() = selection.higgs->GetMomentum().M();
        syncTree.pt_tt()  = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()).Pt();
//        syncTree->m_sv_vegas() = selection.svfitResults.fit_vegas.has_valid_mass
//                ? selection.svfitResults.fit_vegas.mass : default_value;
        syncTree.m_sv() = selection.svfitResult.has_valid_mass
                ? selection.svfitResult.mass : Run2::DefaultFillValueForSyncTree();
        syncTree.pt_sv() = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Pt() : Run2::DefaultFillValueForSyncTree();
        syncTree.eta_sv() = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Eta() : Run2::DefaultFillValueForSyncTree();
        syncTree.phi_sv() = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Phi() : Run2::DefaultFillValueForSyncTree();

        // Kinematic fit
//        syncTree->kinfit_bb_tt_mass() = selection.kinfitResults.mass;
//        syncTree->kinfit_bb_tt_convergence() = selection.kinfitResults.convergence;
//        syncTree->kinfit_bb_tt_chi2() = selection.kinfitResults.chi2;
//        syncTree->kinfit_bb_tt_pull_balance() = selection.kinfitResults.pull_balance;



        // MET
//        const TLorentzVector MET_momentum = MakeLorentzVectorPtEtaPhiM(selection.MET_with_recoil_corrections.pt, 0,
//                                                                       selection.MET_with_recoil_corrections.phi, 0);
//        syncTree->pt_tt_MET() = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()
//                                 + MET_momentum).Pt();

        syncTree.met() = selection.pfMET->Pt();
        syncTree.metphi() = selection.pfMET->Phi();
        syncTree.isPFMET() = selection.pfMET->isPFMET();
        syncTree.metcov00() = selection.pfMET->GetCovVector().at(0);
        syncTree.metcov01() = selection.pfMET->GetCovVector().at(1);
        syncTree.metcov10() = selection.pfMET->GetCovVector().at(2);
        syncTree.metcov11() = selection.pfMET->GetCovVector().at(3);
//        syncTree->mvamet() = MET_momentum.Pt();
//        syncTree->mvametphi() = MET_momentum.Phi();
//        //syncTree->pzetavis();
//        //syncTree->pzetamiss();
//        if(selection.pfMET.significanceMatrix.size()) {
//            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(selection.pfMET.significanceMatrix);
//
//        }
//        const TMatrixD metMVAcov =
//                ntuple::VectorToSignificanceMatrix(selection.MET_with_recoil_corrections.significanceMatrix);
//        syncTree->mvacov00() = metMVAcov[0][0];
//        syncTree->mvacov01() = metMVAcov[0][1];
//        syncTree->mvacov10() = metMVAcov[1][0];
//        syncTree->mvacov11() = metMVAcov[1][1];

        // Leg 1, lepton
        const pat::Muon& patMuon = selection.GetLeg(1)->GetNtupleObject<pat::Muon>();

        syncTree.pt_1()     = selection.GetLeg(1)->GetMomentum().Pt();
        syncTree.phi_1()    = selection.GetLeg(1)->GetMomentum().Phi();
        syncTree.eta_1()    = selection.GetLeg(1)->GetMomentum().Eta();
        syncTree.m_1()      = selection.GetLeg(1)->GetMomentum().M();
        syncTree.q_1()      = selection.GetLeg(1)->GetCharge();
        syncTree.pfmt_1()     = Calculate_MT(selection.GetLeg(1)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree.d0_1()     = Calculate_dxy(selection.GetLeg(1)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(1)->GetMomentum());
        syncTree.dZ_1()     = selection.GetLeg(1)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
        syncTree.iso_1()    = (patMuon.pfIsolationR03().sumChargedHadronPt + std::max(
                                   patMuon.pfIsolationR03().sumNeutralHadronEt +
                                   patMuon.pfIsolationR03().sumPhotonEt -
                                   0.5 * patMuon.pfIsolationR03().sumPUPt, 0.0)) / patMuon.pt();

        syncTree.id_e_mva_nt_loose_1() = Run2::DefaultFillValueForSyncTree();

        syncTree.gen_match_1() = true;
        // Leg 2, tau
        const pat::Tau& patTau = selection.GetLeg(2)->GetNtupleObject<pat::Tau>();

        syncTree.pt_2()     = selection.GetLeg(2)->GetMomentum().Pt();
        syncTree.phi_2()    = selection.GetLeg(2)->GetMomentum().Phi();
        syncTree.eta_2()    = selection.GetLeg(2)->GetMomentum().Eta();
        syncTree.m_2()      = selection.GetLeg(2)->GetMomentum().M();
        syncTree.q_2()      = selection.GetLeg(2)->GetCharge();
        syncTree.pfmt_2()     = Calculate_MT(selection.GetLeg(2)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree.d0_2()     = Calculate_dxy(selection.GetLeg(2)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(2)->GetMomentum());
        syncTree.dZ_2()     = selection.GetLeg(2)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
	syncTree.iso_2()    = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        syncTree.id_e_mva_nt_loose_1() = Run2::DefaultFillValueForSyncTree();
        syncTree.gen_match_2() = true;

        syncTree.againstElectronLooseMVA6_2()   = patTau.tauID("againstElectronLooseMVA6");
        syncTree.againstElectronMediumMVA6_2()  = patTau.tauID("againstElectronMediumMVA6");
        syncTree.againstElectronTightMVA6_2()   = patTau.tauID("againstElectronTightMVA6");
        syncTree.againstElectronVLooseMVA6_2()  = patTau.tauID("againstElectronVLooseMVA6");
        syncTree.againstElectronVTightMVA6_2()  = patTau.tauID("againstElectronVTightMVA6");

        syncTree.againstMuonLoose3_2()          = patTau.tauID("againstMuonLoose3");
        syncTree.againstMuonTight3_2()          = patTau.tauID("againstMuonTight3");

        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = patTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        syncTree.byIsolationMVA3newDMwLTraw_2()               = patTau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        syncTree.byIsolationMVA3oldDMwLTraw_2()               = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        syncTree.byIsolationMVA3newDMwoLTraw_2()              = Run2::DefaultFillValueForSyncTree();
        syncTree.byIsolationMVA3oldDMwoLTraw_2()              = Run2::DefaultFillValueForSyncTree();

        syncTree.decayModeFindingOldDMs_2() = patTau.tauID("decayModeFinding");

        syncTree.dilepton_veto() = selection.Zveto;

        // Jets
        syncTree.njetspt20() = selection.jets.size();
        //syncTree.njets()     = selection.numJet;
        syncTree.nbtag()     = selection.bjets.size();
        //int jetCount = 0;

        Int_t numJet = 0;
        for( const CandidateV2Ptr& jet : selection.jets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                if (jet->GetMomentum().Pt() > 30 ) numJet++;
                syncTree.pt_jets()   .push_back(jet->GetMomentum().Pt());
                syncTree.eta_jets()  .push_back(jet->GetMomentum().Eta());
                syncTree.phi_jets()  .push_back(jet->GetMomentum().Phi());
                syncTree.energy_jets() .push_back(jet->GetMomentum().E());
                syncTree.rawf_jets() .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree.mva_jets()  .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree.csv_jets()  .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
         }
        syncTree.njets() = numJet;


        for( const CandidateV2Ptr& jet : selection.bjets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                syncTree.pt_bjets()   .push_back(jet->GetMomentum().Pt());
                syncTree.eta_bjets()  .push_back(jet->GetMomentum().Eta());
                syncTree.phi_bjets()  .push_back(jet->GetMomentum().Phi());
                syncTree.energy_bjets() .push_back(jet->GetMomentum().E());
                syncTree.rawf_bjets() .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree.mva_bjets()  .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree.csv_bjets()  .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
         }

//        for (const CandidatePtr& jet : selection.bjets_all) {
//            const ntuple::Jet& ntuple_jet = jet->GetNtupleObject<ntuple::Jet>();

//            syncTree->pt_Bjets()      .push_back( jet->GetMomentum().Pt() );
//            syncTree->eta_Bjets()     .push_back( jet->GetMomentum().Eta() );
//            syncTree->phi_Bjets()     .push_back( jet->GetMomentum().Phi() );
//            syncTree->energy_Bjets()  .push_back( jet->GetMomentum().E() );
//            syncTree->chargedHadronEF_Bjets().push_back( ntuple_jet.chargedHadronEnergyFraction );
//            syncTree->neutralHadronEF_Bjets()  .push_back( ntuple_jet.neutralHadronEnergyFraction );
//            syncTree->photonEF_Bjets()         .push_back( ntuple_jet.photonEnergyFraction );
//            syncTree->muonEF_Bjets()  .push_back( ntuple_jet.muonEnergyFraction );
//            syncTree->electronEF_Bjets()  .push_back( ntuple_jet.electronEnergyFraction );
//            syncTree->csv_Bjets()     .push_back( ntuple_jet.combinedSecondaryVertexBJetTags );
//            // inspect the flavour of the gen jet
//            const VisibleGenObjectVector matched_bjets_MC = FindMatchedObjects(jet->GetMomentum(),
//                                                                               selection.GetFinalStateMC().b_jets,
//                                                                               cuts::DeltaR_MC_Match);
//            const bool isJet_MC_Bjet = matched_bjets_MC.size() != 0;
//            const bool isJet_MC_Bjet_withLeptonicDecay = isJet_MC_Bjet
//                    && matched_bjets_MC.at(0).finalStateChargedLeptons.size() != 0;
//            syncTree->isBjet_MC_Bjet()                  .push_back( isJet_MC_Bjet );
//            syncTree->isBjet_MC_Bjet_withLeptonicDecay().push_back( isJet_MC_Bjet_withLeptonicDecay );
//        }

//        syncTree->x_PV() = primaryVertex->GetPosition().x();
//        syncTree->y_PV() = primaryVertex->GetPosition().y();
//        syncTree->z_PV() = primaryVertex->GetPosition().z();
        syncTree.Fill();
    }

//define this as a plug-in
DEFINE_FWK_MODULE(SyncTreeProducer);
