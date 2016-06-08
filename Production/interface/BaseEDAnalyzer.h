/*! Definition of BaseEDAnalyzer class which is the base class for all X->HH->bbTauTau and H->tautau analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <iomanip>
#include <functional>
#include <string>
#include <iostream>


//For CMSSW
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// import LHEEventProduction definition
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//HHbbTauTau Framework
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/CutTools.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/Candidate.h"
#include "h-tautau/Analysis/include/CandidateUtilities.h"
#include "h-tautau/Analysis/include/SyncTree.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"

//SVFit
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

//HHKinFit
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"


#define SELECTION_ENTRY(name) \
    ANA_DATA_ENTRY(cuts::ObjectSelector, name) \
    /**/

#define X(name) \
    object.name

#define Y(name) \
    name

using namespace edm;

namespace analysis {

namespace detail {
const std::map<analysis::Channel,analysis::CandidateV2::Type> channelToTypeMap = {{analysis::Channel::MuTau, analysis::CandidateV2::Type::Muon},
                                                                            {analysis::Channel::ETau, analysis::CandidateV2::Type::Electron},
                                                                            {analysis::Channel::TauTau, analysis::CandidateV2::Type::Tau}};

}

class BaseEDAnalyzerData : public root_ext::AnalyzerData {
public:
    BaseEDAnalyzerData(std::shared_ptr<TFile> outputFile) : AnalyzerData(outputFile) {}
    BaseEDAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};

namespace sv_fit{
struct FitResults {
    static constexpr double default_value = std::numeric_limits<double>::lowest();

    bool has_valid_mass;
    double mass;

    bool has_valid_momentum;
    TLorentzVector momentum;

    FitResults() : has_valid_mass(false), mass(default_value), has_valid_momentum(false) {}
}; //fix SVFit.h
}//sv_fit namespace

struct KinFitResult {
    double mass, chi2, probability;
    int convergence;
    bool HasMass() const { return convergence > 0; }

    KinFitResult() : convergence(std::numeric_limits<int>::lowest()) {}
};

// exported from other .h file
class SelectionManager {
public:
    SelectionManager(root_ext::AnalyzerData& _anaData, const std::string& _selection_label, double _weight)
        : anaData(&_anaData), selection_label(_selection_label), weight(_weight) {}

    template<typename ValueType>
    ValueType FillHistogram(ValueType value, const std::string& histogram_name)
    {
        auto& histogram = anaData->Get(static_cast<TH1D*>(nullptr), histogram_name, selection_label);
        return cuts::fill_histogram(value, histogram, weight);
    }

private:
    root_ext::AnalyzerData* anaData;
    std::string selection_label;
    double weight;
};

struct SelectionResultsV2 {
    static constexpr size_t NumberOfLegs = 2;

    virtual ~SelectionResultsV2() {}
    CandidateV2Ptr higgs;
    Float_t numtruepileupinteractions =-1;
    Float_t HT;
    Int_t NOutPartons;
    Double_t weightevt;
    bool Zveto, electronVeto, muonVeto;
    sv_fit::FitResults svfitResult;
    std::shared_ptr<KinFitResult> kinFit_result;
    //kinematic_fit::four_body::FitResults kinfitResults;
    CandidateV2PtrVector jets;
    int numJet=-1, npv=-1;
    CandidateV2PtrVector jetsTight;
    CandidateV2PtrVector bjets;
    VertexV2PtrVector vertices;
    MissingETPtr pfMET;
    analysis::EventEnergyScale eventEnergyScale;
    //ntuple::MET MET_with_recoil_corrections;
    //ntuple::EventType eventType;

    VertexV2Ptr GetPrimaryVertex() const { return vertices.front(); }
    virtual CandidateV2Ptr GetLeg(size_t leg_id) const = 0;
    //virtual const finalState::bbTauTau& GetFinalStateMC() const = 0;
};
//end export


inline int genMatch( const reco::Candidate::LorentzVector& p4, const std::vector<reco::GenParticle>& genParticles){
  int match = 6;
  float dr = 0.2;
  float dr_new = 0.2;

  for ( std::vector<reco::GenParticle>::const_iterator genp = genParticles.begin();
	genp != genParticles.end(); ++genp ) 
	{
    dr_new = ROOT::Math::VectorUtil::DeltaR( p4, genp->p4());
    
    if (dr_new < dr){
      if (TMath::Abs(genp->pdgId()) == 11) {
	      if (genp->pt() > 8) {
	         if (genp->statusFlags().isPrompt()) {
	            match = 1;
	            dr = dr_new;
	         }
	         else if (genp->statusFlags().isDirectPromptTauDecayProduct()) {
	            match = 3;
	            dr = dr_new;
	         }
	      }
      }
      else if (TMath::Abs(genp->pdgId()) == 13) {
	      if (genp->pt() > 8) {
	         if (genp->isPromptFinalState()){
	            match = 2;
	            dr = dr_new;
	         }
	         if (genp->isDirectPromptTauDecayProductFinalState()){
	            match = 4;
	            dr = dr_new;
	         }
	      }
      }
      else if (TMath::Abs(genp->pdgId()) == 15) {
	      if (genp->statusFlags().isPrompt()) {
	         reco::Candidate::LorentzVector tau_p4 = utils_genMatch::getVisMomentumNoLep(&(*genp));
         	if (tau_p4.pt() > 15) {
	             match = 5;
	            dr = dr_new;
	         }
	      }
      }
    }
  } //GenParticle Loop

  return match;
}



} //namespace analysis


class BaseEDAnalyzer : public edm::EDAnalyzer {

// typedef std::pair<analysis::CandidateV2Ptr&, TLorentzVector&> l1JetMatchPair;
// typedef std::vector<l1JetMatchPair> l1JetMatchPairVtr;

protected:
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
private:
	  edm::EDGetToken electronsMiniAODToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_,eleMediumIdMapToken_,eleCutBasedVetoMapToken_;
    edm::EDGetToken tausMiniAODToken_;
    edm::EDGetToken muonsMiniAODToken_;
    edm::EDGetToken vtxMiniAODToken_;
    edm::EDGetToken pfMETAODToken_;
    edm::EDGetToken jetsMiniAODToken_;
    edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrixTAG_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUInfo_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_;
    edm::EDGetTokenT<GenEventInfoProduct> genWeights_;
    edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT< std::vector<l1extra::L1JetParticle> > l1JetParticle_;
    bool isMC_;
    std::string sampleType;


protected:
  edm::Handle<std::vector<pat::Electron> > electrons;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions, medium_id_decisions, ele_cutBased_veto;
  edm::Handle<std::vector<pat::Tau> > taus;
  edm::Handle<std::vector<pat::Muon> > muons;
  edm::Handle<edm::View<reco::Vertex> > vertices;
  edm::Handle<edm::View<pat::MET> > pfMETs;
  edm::Handle<std::vector<pat::Jet> > jets;
  edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> metCovMatrix;
  edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<LHEEventProduct> lheEventProduct;
  edm::Handle<GenEventInfoProduct> genEvt;
  edm::Handle<std::vector<reco::GenParticle> > pruned;
  edm::Handle< std::vector<l1extra::L1JetParticle> > l1JetParticle;
  


  const bool	isMC() const {return isMC_;}
  const edm::Handle<edm::View<reco::Vertex> > GetVertexCollection() const {return vertices;}
  const edm::Handle<edm::ValueMap<bool> >     GetTElectronIDValueMap() const {return tight_id_decisions;}
  const edm::Handle<edm::ValueMap<bool> >     GetMElectronIDValueMap() const {return medium_id_decisions;}
  const edm::Handle<edm::View<pat::MET> >     GetPFMet()			const {return pfMETs;}
  const edm::Handle<edm::TriggerResults>      GetTriggerBits()      const {return triggerBits;}
  const edm::Handle<GenEventInfoProduct>      GetGenEventInfo()  const {return genEvt;}
  const edm::Handle<std::vector<reco::GenParticle> > GetGenParticles() const {return pruned;}
  const edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> GetMETCovMatrix() {return metCovMatrix;}
  const std::string							  GetSampleType()		const {return sampleType;}

  // static std::shared_ptr<ntuple::SyncTree> GetSyncTree(TFile* file)
//   {
//       static std::shared_ptr<ntuple::SyncTree> sync_tree;
//       if(!sync_tree)
//           sync_tree = std::shared_ptr<ntuple::SyncTree>(new ntuple::SyncTree(file, false));
//       return sync_tree;
//   }

public:
    BaseEDAnalyzer(const edm::ParameterSet& iConfig):
        outputFile(root_ext::CreateRootFile("cuts.root")),
        anaDataBeforeCut(outputFile, "before_cut"), anaDataAfterCut(outputFile, "after_cut"),
        anaDataFinalSelection(outputFile, "final_selection"),
        electronsMiniAODToken_(mayConsume<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  		  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
        eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
        eleCutBasedVetoMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleCutBasedVeto"))),
        tausMiniAODToken_(mayConsume<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
        muonsMiniAODToken_(mayConsume<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
        vtxMiniAODToken_(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
        pfMETAODToken_(mayConsume<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
        jetsMiniAODToken_(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
        metCovMatrixTAG_(consumes<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >>(iConfig.getParameter<edm::InputTag>("metCov"))),
        PUInfo_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
        triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
        triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
        triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
        lheEventProduct_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
        genWeights_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
        prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
        l1JetParticle_(mayConsume< std::vector<l1extra::L1JetParticle> >(iConfig.getParameter<edm::InputTag>("l1JetParticleProduct"))),
        isMC_(iConfig.getParameter<bool>("isMC")),
        sampleType(iConfig.getParameter<std::string>("sampleType")){
     }

    virtual ~BaseEDAnalyzer() {}

    virtual void beginJob() = 0;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) = 0;
    virtual void endJob() = 0;
    //virtual void ProcessEvent(const edm::Event&, const edm::EventSetup&) = 0;

    //  https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    //  PFJetID is tuned on Uncorrected Jet values
    static bool passPFLooseId(const pat::Jet& pat_jet/*, const analysis::CandidateV2Ptr& jet*/)
    {
      //TLorentzVector momentum = jet->GetMomentum();
      const pat::Jet& patJet = pat_jet.correctedJet("Uncorrected");
        //momentum.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
      if(std::abs(patJet.eta())<3.0){
         //if(momentum.E() == 0)                                  return false;
         if(patJet.neutralHadronEnergyFraction() > 0.99)   return false;
         if(patJet.neutralEmEnergyFraction()     > 0.99)   return false;
         if(patJet.nConstituents() <  1)                   return false;
         if(patJet.chargedHadronEnergyFraction() <= 0 && std::abs(patJet.eta()) < 2.4 ) return false;
         if(patJet.chargedEmEnergyFraction() >  0.99  && std::abs(patJet.eta()) < 2.4 ) return false;
         if(patJet.chargedMultiplicity()     <= 0      && std::abs(patJet.eta()) < 2.4 ) return false;
      }
      if(std::abs(patJet.eta())>3.0){
         if(patJet.neutralEmEnergyFraction()     > 0.90)   return false;
         if(patJet.neutralMultiplicity() < 10 )            return false;
      }
      return true;
    }


    void Initialize(const edm::Event& iEvent){
      iEvent.getByToken(electronsMiniAODToken_, electrons);
      iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
      iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
      iEvent.getByToken(eleCutBasedVetoMapToken_,ele_cutBased_veto);
      iEvent.getByToken(tausMiniAODToken_, taus);
      iEvent.getByToken(muonsMiniAODToken_, muons);
      iEvent.getByToken(vtxMiniAODToken_, vertices);
      iEvent.getByToken(pfMETAODToken_, pfMETs);
      iEvent.getByToken(jetsMiniAODToken_, jets);
      iEvent.getByToken(metCovMatrixTAG_,metCovMatrix);
      iEvent.getByToken(triggerBits_, triggerBits);
      iEvent.getByToken(PUInfo_, PUInfo);
      iEvent.getByToken(triggerPrescales_, triggerPrescales);
      iEvent.getByToken(triggerObjects_, triggerObjects);
      iEvent.getByToken(l1JetParticle_,l1JetParticle);
      if (isMC()) iEvent.getByToken(genWeights_,genEvt);
      iEvent.getByToken(prunedGenToken_,pruned);
      try {iEvent.getByToken(lheEventProduct_, lheEventProduct);} catch(...){;};
    }


protected:

    virtual analysis::BaseEDAnalyzerData& GetAnaData() = 0;
    virtual analysis::CandidateV2Ptr SelectHiggs(analysis::CandidateV2PtrVector& higgses) = 0;
    virtual analysis::Channel ChannelId() const = 0;

    analysis::CandidateV2PtrVector ApplyL1TriggerTauMatch(const analysis::CandidateV2PtrVector& higgses){
    
    	analysis::CandidateV2PtrVector L1Higgses;
    	for (const auto& higgs : higgses){
//    	    std::cout<<"***********Higgs**********"<<std::endl;
    		if (L1TauMatch(*higgs)) {
//    			std::cout<<"  Higgs Chosen =  "<<higgs->GetMomentum()<<std::endl;
    			L1Higgses.push_back(higgs);}
    	}
    	return L1Higgses;
    }
    
    
    bool L1TauMatch(const analysis::CandidateV2& candidate){
    
    	std::vector<size_t> pairVector;
    	
    		if(candidate.GetFinalStateDaughters().size()) {
            	for(const auto& daughter : candidate.GetFinalStateDaughters()) {
//            		std::cout<<"------------------------"<<std::endl;
            		size_t jetPos = 0;
            		for(const auto& isoTau : *l1JetParticle){
            			jetPos+=1;
            			TLorentzVector l1JetMomentum;
            			l1JetMomentum.SetPtEtaPhiE(isoTau.pt(), isoTau.eta(), isoTau.phi(), isoTau.energy());
            			if (!(isoTau.pt() > 28)) continue;
//            			std::cout<<"Jet position = "<<jetPos<<"  DeltaR = "<<l1JetMomentum.DeltaR(daughter->GetMomentum())<<std::endl;
            			if(l1JetMomentum.DeltaR(daughter->GetMomentum()) < 0.5 ) pairVector.push_back(jetPos);
            		}
    			}
    			// std::cout<<"\t Size vettore  :  "<<pairVector.size()<<std::endl;
//     			const size_t startingSize = pairVector.size();
//     			std::sort (pairVector.begin(),pairVector.end());
//     			std::vector<size_t>::iterator it;
//     			it = std::unique (pairVector.begin(),pairVector.end());
//     			pairVector.erase(it,pairVector.end());
//     			std::cout<<"\t Size vettore dopo unique :  "<<pairVector.size()<< "   prima:  "<<startingSize<<std::endl;
//     			if (startingSize > 1 && pairVector.size() == startingSize  ) return true;
    			if (pairVector.size() == 2 && pairVector.at(0)!=pairVector.at(1)) return true;
    		}
    		
    	
    	return false;
    }
    
    bool HaveTriggerFired(const edm::Event& iEvent,const std::set<std::string>& hltPaths){

    	  const edm::TriggerNames &triggerNames= iEvent.triggerNames(*(GetTriggerBits()));
//          std::cout << "\n === TRIGGER PATHS === " << std::endl;

          for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
              for (const std::string& triggerPath : hltPaths ){

                  const std::string& objectMatchedPath = triggerNames.triggerName(i);
                  size_t found = objectMatchedPath.find(triggerPath);


                  if(found != std::string::npos){
//                     std::cout << "TriggerPath --->   " << triggerPath << std::endl;
//                     std::cout << "Trigger " << triggerNames.triggerName(i) <<
//                     ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
//                     ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
//                     << std::endl;
                      if(triggerPrescales->getPrescaleForIndex(i) == 1 && triggerBits->accept(i)) return true;
                }
              }
          }
          return false;
    }


    analysis::CandidateV2PtrVector FindCompatibleObjects(const analysis::CandidateV2PtrVector& objects1,
                                                         const analysis::CandidateV2PtrVector& objects2,
                                                         double minDeltaR, analysis::CandidateV2::Type type,
                                                         const std::string& hist_name, int expectedCharge=analysis::CandidateV2::UnknownCharge())
        {
            analysis::CandidateV2PtrVector result;
            for(const analysis::CandidateV2Ptr& object1 : objects1) {
                for(const analysis::CandidateV2Ptr& object2 : objects2) {
                    // std::cout << " \t\t\t  Pair DeltaR = "
                    //           << object2->GetMomentum().DeltaR(object1->GetMomentum()) << std::endl;
                    if(object2->GetMomentum().DeltaR(object1->GetMomentum()) > minDeltaR) {
                        const analysis::CandidateV2Ptr candidate(new analysis::CandidateV2(type, object1, object2));
                        if (expectedCharge != analysis::CandidateV2::UnknownCharge()
                                && candidate->GetCharge() != expectedCharge) continue;
                        result.push_back(candidate);
                        GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),1);
                    }
                }
            }

            GetAnaData().N_objects(hist_name).Fill(result.size(),1);
            return result;
        }

    analysis::CandidateV2PtrVector ApplyTriggerMatch(const edm::Event& iEvent,
                                                     const analysis::CandidateV2PtrVector& higgses,
                                                     const std::set<std::string>& hltPaths,
                                                     bool useStandardTriggerMatch, const bool isCrossTrigger)
    {
        analysis::CandidateV2PtrVector triggeredHiggses;
    	const edm::TriggerNames &triggerNames= iEvent.triggerNames(*triggerBits);

        for (const auto& higgs : higgses){

//            std::cout<<"### Higgs Pair :  \n"<<higgs->GetMomentum()<<std::endl;
            if(!useStandardTriggerMatch && HaveTriggerMatched(*triggerObjects,triggerNames, hltPaths, *higgs,
                                                                        cuts::Htautau_2015::DeltaR_triggerMatch, isCrossTrigger)){
//                std::cout<<"### Pushed Back "<<std::endl;
                triggeredHiggses.push_back(higgs);
            }
    //        if (useStandardTriggerMatch && SyncTreeProducer::HaveTriggerMatched(hltPaths, *higgs))
    //            triggeredHiggses.push_back(higgs);
        }
        return triggeredHiggses;
    }

    bool HaveTriggerMatched(pat::TriggerObjectStandAloneCollection triggerObjects, const edm::TriggerNames &names,
                            const std::set<std::string>& interestingPaths, const analysis::CandidateV2& candidate,
                            const double deltaR_Limit, const bool isCrossTrigger)
    {
        for (const std::string& interestinPath : interestingPaths){
//            std::cout<<"Trigger Path  :   "<<interestinPath<<std::endl;
            if (BaseEDAnalyzer::HaveTriggerMatched(triggerObjects,names,interestinPath,candidate, deltaR_Limit, isCrossTrigger)) return true;
        }
        return false;
    }
//Now it could manage also the single object Trigger
    bool HaveTriggerMatched(pat::TriggerObjectStandAloneCollection triggerObjects, const edm::TriggerNames &names,
                            const std::string& interestingPath, const analysis::CandidateV2& candidate,
                            const double deltaR_Limit, const bool isCrossTrigger)
    {
        if(candidate.GetFinalStateDaughters().size()) {
            for(const auto& daughter : candidate.GetFinalStateDaughters()) {
                if(isCrossTrigger &&
                        !BaseEDAnalyzer::HaveTriggerMatched(triggerObjects, names, interestingPath, *daughter, deltaR_Limit,isCrossTrigger))
                    return false;
                if(!isCrossTrigger &&
                        daughter->GetType() == analysis::detail::channelToTypeMap.at(ChannelId()) &&
                        BaseEDAnalyzer::HaveTriggerMatched(triggerObjects, names, interestingPath, *daughter, deltaR_Limit,isCrossTrigger))
                    return true;
            }
            if ( isCrossTrigger)return true;
            if ( !isCrossTrigger)return false;
        }

        for (auto &triggerObject : triggerObjects){
            triggerObject.unpackPathNames(names);
            TLorentzVector triggerObjectMomentum;
            triggerObjectMomentum.SetPtEtaPhiM(triggerObject.pt(), triggerObject.eta(), triggerObject.phi(), triggerObject.mass());

            for (unsigned n = 0; n < triggerObject.pathNames(true).size(); ++n){
                const std::string& objectMatchedPath = triggerObject.pathNames(true).at(n);

                const size_t found = objectMatchedPath.find(interestingPath);
                bool isBoth = triggerObject.hasPathName( triggerObject.pathNames(true).at(n), true, true );
                if (found != std::string::npos) {
//                    std::cout<<"\t Trigger obj :   "<<triggerObjectMomentum<<std::endl;
//                    std::cout << "\t\t Path Names :   " << objectMatchedPath
//                          << (isBoth ? "  PASS\n " : "  fail (or not run) --\n ")
//                          //<< triggerObject.pathLastFilterAccepted.at(n) << "\n"
//                          << "\t\t Candidato :  " << candidate.GetMomentum()
//                          <<"\t\t DeltaR = "<< triggerObjectMomentum.DeltaR(candidate.GetMomentum())
//                          <<(triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit ? "\t Matched":"\t NotMatched")
//                          << std::endl;
                }
                if (found != std::string::npos && isBoth &&
                        triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit)
                    return true;
            }
        }
        return false;
    }

    std::pair<double,int> computeHtValue (){
        if(lheEventProduct.isValid()){
            const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
            std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
            double lheHt = 0.;
            int lheNOutPartons = 0;
            size_t numParticles = lheParticles.size();
            for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
                int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
                int status = lheEvent.ISTUP[idxParticle];
                if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
                    lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
                    ++lheNOutPartons;
                }
            }
            std::pair<double,int> lheInfoPair(lheHt,lheNOutPartons);
            return lheInfoPair;
        }
        else {
          std::pair<double,int> lheInfoPair(-999.99,-1);
          return lheInfoPair;
        }
    }


    template<typename FirstLegPATType>
    analysis::sv_fit::FitResults SVFit(const analysis::CandidateV2Ptr& higgs, const analysis::MissingETPtr& met){

        static const std::map<analysis::CandidateV2::Type, svFitStandalone::kDecayType> decayTypeMap = {
            { analysis::CandidateV2::Type::Electron, svFitStandalone::kTauToElecDecay },
            { analysis::CandidateV2::Type::Muon, svFitStandalone::kTauToMuDecay },
            { analysis::CandidateV2::Type::Tau, svFitStandalone::kTauToHadDecay }
        };

        //const FirstLegPATType& leg1 = higgs->GetDaughters().at(0)->GetNtupleObject<FirstLegPATType>();
        const TLorentzVector leg1_momentum = higgs->GetDaughters().at(0)->GetMomentum();
        const pat::Tau&  tau  = higgs->GetDaughters().at(1)->GetNtupleObject<pat::Tau>();
        const TLorentzVector tau_momentum = higgs->GetDaughters().at(1)->GetMomentum();
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
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(decayTypeMap.at(higgs->GetDaughters().at(0)->GetType()), leg1_momentum.Pt(),
                                                                          leg1_momentum.Eta(), leg1_momentum.Phi(), leg1_momentum.M())); // tau -> electron decay (Pt, eta, phi, mass)
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(decayTypeMap.at(higgs->GetDaughters().at(1)->GetType()), tau_momentum.Pt(),
                                                                          tau_momentum.Eta(), tau_momentum.Phi(), tau_momentum.M(),tau.decayMode())); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
          SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
          algo.addLogM(false);
          edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
          TH1::AddDirectory(false);
          TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
          algo.shiftVisPt(true, inputFile_visPtResolution);
          algo.integrateMarkovChain();

          analysis::sv_fit::FitResults svfitResult;
          //double mass = algo.getMass(); // return value is in units of GeV
          if ( algo.isValidSolution() ) {
              svfitResult.mass = algo.mass();
              svfitResult.has_valid_mass = true;
              //if(fitAlgorithm == FitAlgorithm::MarkovChain) {
                  svfitResult.momentum.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.mass());
                  svfitResult.has_valid_momentum = true;
             // }
//             std::cout << "... m svfit : " << algo.mass() << " +/- " << algo.massUncert() << std::endl;
          } else {
//            std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
          }

          delete inputFile_visPtResolution;

          return svfitResult;
    }

   void FillSyncTree(const edm::Event& iEvent, const analysis::SelectionResultsV2& selection, ntuple::SyncTree& syncTree)
       {
           using namespace analysis;
        static const float default_value = Run2::DefaultFloatFillValueForSyncTree();


		    const VertexV2Ptr primaryVertex(new VertexV2((*(BaseEDAnalyzer::GetVertexCollection())).front()));

        //BaseEDAnalyzer::FillSyncTree(iEvent, selection, syncTree);

        //const auto primaryVertex = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);

        // Event
//        std::cout<<"~~~~~~~~~~~~~EVENT Info~~~~~~~~~"<<std::endl;
//        std::cout<<"Run = "<<iEvent.id().run()<<"  Lumi = "<<iEvent.id().luminosityBlock()
//                <<" Event = "<<iEvent.id().event()<<std::endl;
        syncTree().run  = iEvent.id().run();
        syncTree().lumi = iEvent.id().luminosityBlock();
        syncTree().evt  = iEvent.id().event();
        syncTree().channelID = static_cast<int>(ChannelId());
        syncTree().eventEnergyScale = static_cast<int>(selection.eventEnergyScale);


        syncTree().HT   = selection.HT;
        syncTree().NOutPartons = selection.NOutPartons;
        if(selection.HT<100) syncTree().HTBin = static_cast<int>(Run2::HTbinning::lt100);
        if(100<=selection.HT && selection.HT<200) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f100to200);
        if(200<=selection.HT && selection.HT<400) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f200to400);
        if(400<=selection.HT && selection.HT<600) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f400to600);
        if(selection.HT>=600) syncTree().HTBin = static_cast<int>(Run2::HTbinning::gt600);
        if(selection.HT<0)    syncTree().HTBin = -1;

         if (BaseEDAnalyzer::isMC()) syncTree().weightevt = selection.weightevt;
         else syncTree().weightevt = default_value;

        syncTree().npv = selection.npv;
        syncTree().npu = selection.numtruepileupinteractions;


        // HTT candidate
        syncTree().m_vis = selection.higgs->GetMomentum().M();
        syncTree().pt_tt  = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()).Pt();
        syncTree().m_sv = selection.svfitResult.has_valid_mass
                ? selection.svfitResult.mass : Run2::DefaultFillValueForSyncTree();
        syncTree().pt_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Pt() : Run2::DefaultFillValueForSyncTree();
        syncTree().eta_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Eta() : Run2::DefaultFillValueForSyncTree();
        syncTree().phi_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Phi() : Run2::DefaultFillValueForSyncTree();
                
        // const pat::MET& patMET = selection.pfMET->GetNtupleObject<pat::MET>();
         
        syncTree().met        = selection.pfMET->Pt();
        syncTree().metphi     = selection.pfMET->Phi();
        syncTree().isPFMET    = selection.pfMET->isPFMET();
        syncTree().metcov00   = selection.pfMET->GetCovVector().at(0);
        syncTree().metcov01   = selection.pfMET->GetCovVector().at(1);
        syncTree().metcov10   = selection.pfMET->GetCovVector().at(2);
        syncTree().metcov11   = selection.pfMET->GetCovVector().at(3);
        
        
        // Jets
        syncTree().njetspt20 = selection.jets.size();
        //syncTree.njets()     = selection.numJet;
        syncTree().nbtag     = selection.bjets.size();
        //int jetCount = 0;

        Int_t numJet = 0;
        for( const CandidateV2Ptr& jet : selection.jets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                if (jet->GetMomentum().Pt() > 30 ) numJet++;
                syncTree().pt_jets     .push_back(jet->GetMomentum().Pt());
                syncTree().eta_jets    .push_back(jet->GetMomentum().Eta());
                syncTree().phi_jets    .push_back(jet->GetMomentum().Phi());
                syncTree().energy_jets .push_back(jet->GetMomentum().E());
                syncTree().rawf_jets   .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree().mva_jets    .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree().csv_jets    .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                //std::cout<<" \t\t CSV:  "<< pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
                syncTree().partonFlavour_jets .push_back(pat_jet1.partonFlavour());
         }
        syncTree().njets = numJet;

        for( const CandidateV2Ptr& jet : selection.bjets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                syncTree().pt_bjets     .push_back(jet->GetMomentum().Pt());
                syncTree().eta_bjets    .push_back(jet->GetMomentum().Eta());
                syncTree().phi_bjets    .push_back(jet->GetMomentum().Phi());
                syncTree().energy_bjets .push_back(jet->GetMomentum().E());
                syncTree().rawf_bjets   .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree().mva_bjets    .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree().csv_bjets    .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                syncTree().partonFlavour_bjets .push_back(pat_jet1.partonFlavour());
         }
         
         // syncTree().kinFit_m = selection.kinFit_result->mass;
//          syncTree().kinFit_m = selection.kinFit_result->chi2;
//          syncTree().kinFit_m = selection.kinFit_result->probability;
//          syncTree().kinFit_m = selection.kinFit_result->convergence;
         
         syncTree().dilepton_veto  = selection.Zveto;
        syncTree().extraelec_veto = selection.electronVeto;
        syncTree().extramuon_veto = selection.muonVeto;
        
   }


    template<typename ObjectType, typename BaseSelectorType, typename NtupleObjectType,
             typename ObjectPtrType = std::shared_ptr<const ObjectType>,
             typename Comparitor = std::less<ObjectPtrType> >
    std::vector<ObjectPtrType> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                              const std::vector<NtupleObjectType>& ntuple_objects,
                                              Comparitor comparitor = Comparitor())
    {
        cuts::ObjectSelector& objectSelector = GetAnaData().Selection(selection_label);
        analysis::SelectionManager selectionManager(anaDataBeforeCut, selection_label, 1);

        const auto selector = [&](size_t id) -> ObjectPtrType {
           // ObjectPtrType candidate(new ObjectType(ntuple_objects.at(id)));
            ObjectPtrType candidate(new ObjectType(ntuple_objects.at(id),id));
            cuts::Cutter cut(&objectSelector);
            base_selector(candidate, selectionManager, cut);
            return candidate;
        };

        const auto selected = objectSelector.collect_objects<ObjectPtrType>(1,
                                                                            ntuple_objects.size(), selector,
                                                                            comparitor);
        // SelectionManager selectionManager_afterCut(anaDataAfterCut, selection_label,
        //                                            GetEventWeights().GetPartialWeight());
        // for(const auto& candidate : selected) {
        //     cuts::Cutter cut(nullptr);
        //     base_selector(candidate, selectionManager_afterCut, cut);
        // }
        // GetAnaData().N_objects(selection_label).Fill(selected.size(), GetEventWeights().GetPartialWeight());
        // GetAnaData().N_objects(selection_label + "_ntuple").Fill(ntuple_objects.size(),
        //                                                          GetEventWeights().GetPartialWeight());
        return selected;
    }

    template<typename BaseSelectorMethod, typename PATObjectType>
    analysis::CandidateV2PtrVector CollectCandidateObjects(const std::string& selection_label, BaseSelectorMethod selector_method,
                                               const std::vector<PATObjectType>& ntuple_objects)
    {
        const auto base_selector = [&](const analysis::CandidateV2Ptr& candidate, analysis::SelectionManager& selectionManager,
                                       cuts::Cutter& cut)
            { (this->*selector_method)(candidate, selectionManager, cut); };
        return CollectObjects<analysis::CandidateV2>(selection_label, base_selector, ntuple_objects);
    }


    analysis::CandidateV2PtrVector CollectSignalMuons()
    {
        return CollectCandidateObjects("muons_sgn", &BaseEDAnalyzer::SelectSignalMuon, *muons);
    }

    virtual void SelectSignalMuon(const analysis::CandidateV2Ptr& muon, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Signal muon selection for signal not implemented");
    }

    analysis::CandidateV2PtrVector CollectSignalElectrons()
    {
        return CollectCandidateObjects("electrons_sgn", &BaseEDAnalyzer::SelectSignalElectron, *electrons);
    }

    virtual void SelectSignalElectron(const analysis::CandidateV2Ptr& electron, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Signal electron selection for signal not implemented");
    }

    analysis::CandidateV2PtrVector CollectZmuons()
    {
        return CollectCandidateObjects("muons", &BaseEDAnalyzer::SelectZMuon, *muons);
    }

    virtual void SelectZMuon(const analysis::CandidateV2Ptr& muon, analysis::SelectionManager& selectionManager,
                             cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015::MuTau;
        const pat::Muon& object = muon->GetNtupleObject<pat::Muon>();
		const auto primaryVertex = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);

		const double iso_mu = (object.pfIsolationR03().sumChargedHadronPt + std::max(
                            object.pfIsolationR03().sumNeutralHadronEt +
                            object.pfIsolationR03().sumPhotonEt -
                            0.5 * object.pfIsolationR03().sumPUPt, 0.0)) / object.pt();

        cut(true, ">0 mu cand");
        cut(X(pt()) > ZmumuVeto::pt, "pt");
        cut(std::abs( X(eta()) ) < ZmumuVeto::eta, "eta");
        const double muonDZ = std::abs(object.muonBestTrack()->dz(primaryVertex->position()));
        cut(Y(muonDZ)  < muonID::dz, "dz");
        const double muonDB = std::abs(object.muonBestTrack()->dxy(primaryVertex->position()));
        cut(Y(muonDB) < muonID::dB, "dxy");
        cut(X(isTrackerMuon()), "trackerMuon");
        cut(X(isGlobalMuon()), "GlobalMuon");
        cut(X(isPFMuon()), "PFMuon");
        cut(Y(iso_mu) < 0.3, "pFRelIso");
    }
    
    analysis::CandidateV2PtrVector CollectZelectrons()
    {
        return CollectCandidateObjects("muons", &BaseEDAnalyzer::SelectZElectron, *electrons);
    }

    virtual void SelectZElectron(const analysis::CandidateV2Ptr& electron, analysis::SelectionManager& selectionManager,
                             cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015::ETau;
        const pat::Electron& object = electron->GetNtupleObject<pat::Electron>();

        const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);
        const Ptr<pat::Electron> elPtr(electrons, electron->GetIndex() );
        const double iso = (object.pfIsolationVariables().sumChargedHadronPt
                           + std::max(object.pfIsolationVariables().sumNeutralHadronEt
                           + object.pfIsolationVariables().sumPhotonEt -
                           0.5 * object.pfIsolationVariables().sumPUPt, 0.0)) / object.pt();


        cut(true, ">0 mu cand");
//    std::cout<<"  edm::Ptr  id  =  "<<elPtr.id()<<"  edm::Ptr  key  =  "<<elPtr.key()<<
//               " MVA TightID  = "<<isTight<<std::endl;

        cut(X(pt()) > ZeeVeto::pt, "pt");
        cut(std::abs( X(eta()) ) < ZeeVeto::eta, "eta");
        const double electronD0 = fabs(object.gsfTrack()->dxy(PV->position()));
        cut(Y(electronD0) < ZeeVeto::d0, "dxy");
        const double electronDZ = fabs(object.gsfTrack()->dz(PV->position()));
        cut(Y(electronDZ) < ZeeVeto::dz, "dz");
        const bool veto  = (*ele_cutBased_veto)[elPtr];
        cut(Y(veto), "cut based veto");
        cut(Y(iso)< ZeeVeto::pfRelIso, "iso");
    }


    analysis::CandidateV2PtrVector CollectSignalTaus()
    {
        return CollectCandidateObjects("taus_sgn", &BaseEDAnalyzer::SelectSignalTau, *taus);
    }

    virtual void SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Signal tau selection for signal not implemented");
    }


    analysis::CandidateV2PtrVector CollectVetoElectrons()
    {
        return CollectCandidateObjects("electrons", &BaseEDAnalyzer::SelectVetoElectron, *electrons);
    }


    analysis::CandidateV2PtrVector CollectVetoMuons()
    {
        return CollectCandidateObjects("muons", &BaseEDAnalyzer::SelectVetoMuon, *muons);
    }

    virtual void SelectVetoMuon(const analysis::CandidateV2Ptr& muon, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015;
        const pat::Muon& object = muon->GetNtupleObject<pat::Muon>();
        const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);
        const double pfIso = (object.pfIsolationR03().sumChargedHadronPt + std::max(
                                  object.pfIsolationR03().sumNeutralHadronEt +
                                  object.pfIsolationR03().sumPhotonEt -
                                  0.5 * object.pfIsolationR03().sumPUPt, 0.0)) / object.pt();

        cut(true, ">0 mu cand");
        cut(X(pt()) > muonVeto::pt, "pt");
        cut(std::abs( X(eta()) ) < muonVeto::eta, "eta");
        const double muonDB = fabs(object.muonBestTrack()->dxy(PV->position()));
        cut(Y(muonDB) < muonVeto::dB, "dxy");
        const double muonDZ = fabs(object.muonBestTrack()->dz(PV->position()));
        cut(Y(muonDZ) < muonVeto::dz, "dz");
        cut(Y(pfIso) < muonVeto::pfRelIso, "iso");
        cut(X(isMediumMuon()), "muonID");
    }

    virtual void SelectVetoElectron(const analysis::CandidateV2Ptr& electron, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015;
        const pat::Electron& object = electron->GetNtupleObject<pat::Electron>();

        const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);
        const Ptr<pat::Electron> elPtr(electrons, electron->GetIndex() );

        double iso = (object.pfIsolationVariables().sumChargedHadronPt
                       + std::max(object.pfIsolationVariables().sumNeutralHadronEt
                       + object.pfIsolationVariables().sumPhotonEt -
                       0.5 * object.pfIsolationVariables().sumPUPt, 0.0)) / object.pt();

        cut(true, ">0 mu cand");
        cut(X(pt()) > electronVeto::pt, "pt");
        cut(std::abs( X(eta()) ) < electronVeto::eta_high, "eta");
        const double electronD0 = fabs(object.gsfTrack()->dxy(PV->position()));
        cut(Y(electronD0) < electronVeto::d0, "dxy");
        const double electronDZ = fabs(object.gsfTrack()->dz(PV->position()));
        cut(Y(electronDZ) < electronVeto::dz, "dz");
        const bool isTight  = (*(BaseEDAnalyzer::GetMElectronIDValueMap()))[elPtr];
        cut(isTight, "electronMVATightID");
        const int eleMissingHits = object.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        cut(eleMissingHits <= electronVeto::missingHits,"missingHits");
        cut(X(passConversionVeto()),"conversionVeto");
        cut(Y(iso) < electronVeto::pFRelIso, "iso");
    }

    analysis::CandidateV2PtrVector CollectJets()
    {
        return CollectCandidateObjects("jets", &BaseEDAnalyzer::SelectJets, *jets);
    }

    virtual void SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
            throw std::runtime_error("Jets selection not implemented");
    }
    
    analysis::CandidateV2PtrVector CollectJetsTight()
    {
        return CollectCandidateObjects("jets", &BaseEDAnalyzer::SelectJetsTight, *jets);
    }

    virtual void SelectJetsTight(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
            throw std::runtime_error("Jets selection not implemented");
    }

    analysis::CandidateV2PtrVector CollectBJets()
    {
        return CollectCandidateObjects("bjets", &BaseEDAnalyzer::SelectBJets, *jets);
    }

    virtual void SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
            throw std::runtime_error("BJets selection not implemented");
    }
};
