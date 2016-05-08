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
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
#include "h-tautau/Analysis/include/SyncTree.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"

//SVFit
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#define SELECTION_ENTRY(name) \
    ANA_DATA_ENTRY(cuts::ObjectSelector, name) \
    /**/

#define X(name) \
    object.name

#define Y(name) \
    name

using namespace edm;

// namespace analysis {
// 	class CandidateV2;
// 	typedef std::shared_ptr<const CandidateV2> CandidateV2Ptr;
// 	typedef std::vector<CandidateV2Ptr> CandidateV2PtrVector;
// 	enum class CandidateV2::Type;
//
// 	class VertexV2;
// 	typedef std::shared_ptr<const VertexV2> VertexV2Ptr;
// 	typedef std::vector<VertexV2Ptr> VertexV2PtrVector;
//
// 	class MissingET;
// 	typedef std::shared_ptr<const MissingET> MissingETPtr;
// };

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
    Double_t weightevt;
    bool Zveto, electronVeto, muonVeto;
    sv_fit::FitResults svfitResult;
    //kinematic_fit::four_body::FitResults kinfitResults;
    CandidateV2PtrVector jets;
    int numJet=-1, npv=-1;
    CandidateV2PtrVector jetsPt20;
    CandidateV2PtrVector bjets;
    CandidateV2PtrVector retagged_bjets;
    VertexV2PtrVector vertices;
    MissingETPtr pfMET;
    //ntuple::MET MET_with_recoil_corrections;
    //ntuple::EventType eventType;

    VertexV2Ptr GetPrimaryVertex() const { return vertices.front(); }
    virtual CandidateV2Ptr GetLeg(size_t leg_id) const = 0;
    //virtual const finalState::bbTauTau& GetFinalStateMC() const = 0;
};
//end export

} //namespace analysis


class BaseEDAnalyzer : public edm::EDAnalyzer {

// typedef std::pair<analysis::CandidateV2Ptr&, TLorentzVector&> l1JetMatchPair;
// typedef std::vector<l1JetMatchPair> l1JetMatchPairVtr;

protected:
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
private:
	  edm::EDGetToken electronsMiniAODToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_,eleMediumIdMapToken_;
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
    edm::EDGetTokenT< std::vector<l1extra::L1JetParticle> > l1JetParticle_;
    bool isMC_;
    bool computeHT_;
    std::string sampleType;

protected:
  edm::Handle<std::vector<pat::Electron> > electrons;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions, medium_id_decisions;
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
  edm::Handle< std::vector<l1extra::L1JetParticle> > l1JetParticle;


  const bool	isMC() const {return isMC_;}
  const edm::Handle<edm::View<reco::Vertex> > GetVertexCollection() const {return vertices;}
  const edm::Handle<edm::ValueMap<bool> >     GetTElectronIDValueMap() const {return tight_id_decisions;}
  const edm::Handle<edm::ValueMap<bool> >     GetMElectronIDValueMap() const {return medium_id_decisions;}
  const edm::Handle<edm::View<pat::MET> >     GetPFMet()			const {return pfMETs;}
  const edm::Handle<edm::TriggerResults>      GetTriggerBits()      const {return triggerBits;}
  const edm::Handle<GenEventInfoProduct>      GetGenEventInfo()  const {return genEvt;}
  const edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >> GetMETCovMatrix() {return metCovMatrix;}
  const std::string							  GetSampleType()		const {return sampleType;}

  static std::shared_ptr<ntuple::SyncTree> GetSyncTree(TFile* file)
  {
      static std::shared_ptr<ntuple::SyncTree> sync_tree;
      if(!sync_tree)
          sync_tree = std::shared_ptr<ntuple::SyncTree>(new ntuple::SyncTree(file, false));
      return sync_tree;
  }

public:
    BaseEDAnalyzer(const edm::ParameterSet& iConfig):
        outputFile(root_ext::CreateRootFile("cuts.root")),
        anaDataBeforeCut(outputFile, "before_cut"), anaDataAfterCut(outputFile, "after_cut"),
        anaDataFinalSelection(outputFile, "final_selection"),
        electronsMiniAODToken_(mayConsume<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  		  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
        eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
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
        l1JetParticle_(mayConsume< std::vector<l1extra::L1JetParticle> >(iConfig.getParameter<edm::InputTag>("l1JetParticleProduct"))),
        isMC_(iConfig.getParameter<bool>("isMC")),
        computeHT_(iConfig.getParameter<bool>("HTBinning")),
        sampleType(iConfig.getParameter<std::string>("sampleType")){
     }

    virtual ~BaseEDAnalyzer() {}

    virtual void beginJob() = 0;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) = 0;
    virtual void endJob() = 0;
    //virtual void ProcessEvent(const edm::Event&, const edm::EventSetup&) = 0;

    //  https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    static bool passPFLooseId(const pat::Jet& jet)
    {
        TLorentzVector momentum;
        momentum.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
        if(momentum.E() == 0)                                  return false;
        if(jet.neutralHadronEnergyFraction() > 0.99)   return false;
        if(jet.neutralEmEnergyFraction()     > 0.99)   return false;
        if(jet.nConstituents() <  1)                   return false;
        if(jet.chargedHadronEnergyFraction() <= 0 && std::abs(jet.eta()) < 2.4 ) return false;
        if(jet.chargedEmEnergyFraction() >  0.99  && std::abs(jet.eta()) < 2.4 ) return false;
        if(jet.chargedMultiplicity()     <= 0      && std::abs(jet.eta()) < 2.4 ) return false;
        return true;
    }


    void Initialize(const edm::Event& iEvent){
      iEvent.getByToken(electronsMiniAODToken_, electrons);
      iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
      iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
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
      try {iEvent.getByToken(lheEventProduct_, lheEventProduct);} catch(...){;};
    }


protected:

    virtual analysis::BaseEDAnalyzerData& GetAnaData() = 0;
    virtual analysis::CandidateV2Ptr SelectHiggs(analysis::CandidateV2PtrVector& higgses) = 0;
    virtual analysis::Channel ChannelId() const = 0;

    analysis::CandidateV2PtrVector ApplyL1TriggerTauMatch(const analysis::CandidateV2PtrVector& higgses){
    
    	analysis::CandidateV2PtrVector L1Higgses;
    	for (const auto& higgs : higgses){
    	    std::cout<<"***********Higgs**********"<<std::endl;
    		if (L1TauMatch(*higgs)) {
    			std::cout<<"  Higgs Chosen =  "<<higgs->GetMomentum()<<std::endl;
    			L1Higgses.push_back(higgs);}
    	}
    	return L1Higgses;
    }
    
    
    bool L1TauMatch(const analysis::CandidateV2& candidate){
    
    	std::vector<size_t> pairVector;
    	
    		if(candidate.GetFinalStateDaughters().size()) {
            	for(const auto& daughter : candidate.GetFinalStateDaughters()) {
            		std::cout<<"------------------------"<<std::endl;
            		size_t jetPos = 0;
            		for(const auto& isoTau : *l1JetParticle){
            			jetPos+=1;
            			TLorentzVector l1JetMomentum;
            			l1JetMomentum.SetPtEtaPhiE(isoTau.pt(), isoTau.eta(), isoTau.phi(), isoTau.energy());
            			if (!(isoTau.pt() > 28)) continue;
            			std::cout<<"Jet position = "<<jetPos<<"  DeltaR = "<<l1JetMomentum.DeltaR(daughter->GetMomentum())<<std::endl;
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
          std::cout << "\n === TRIGGER PATHS === " << std::endl;

          for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
              for (const std::string& triggerPath : hltPaths ){

                  const std::string& objectMatchedPath = triggerNames.triggerName(i);
                  size_t found = objectMatchedPath.find(triggerPath);


                  if(found != std::string::npos){
                     std::cout << "TriggerPath --->   " << triggerPath << std::endl;
                     std::cout << "Trigger " << triggerNames.triggerName(i) <<
                     ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                     ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                     << std::endl;
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

            std::cout<<"### Higgs Pair :  \n"<<higgs->GetMomentum()<<std::endl;
            if(!useStandardTriggerMatch && HaveTriggerMatched(*triggerObjects,triggerNames, hltPaths, *higgs,
                                                                        cuts::Htautau_2015::DeltaR_triggerMatch, isCrossTrigger)){
                std::cout<<"### Pushed Back "<<std::endl;
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
            std::cout<<"Trigger Path  :   "<<interestinPath<<std::endl;
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

        const FirstLegPATType& leg1 = higgs->GetDaughters().at(0)->GetNtupleObject<FirstLegPATType>();
        const pat::Tau&  tau  = higgs->GetDaughters().at(1)->GetNtupleObject<pat::Tau>();
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
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(decayTypeMap.at(higgs->GetDaughters().at(0)->GetType()), leg1.pt(),
                                                                          leg1.eta(), leg1.phi(), leg1.mass())); // tau -> electron decay (Pt, eta, phi, mass)
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(decayTypeMap.at(higgs->GetDaughters().at(1)->GetType()), tau.pt(),
                                                                          tau.eta(), tau.phi(), tau.mass(),tau.decayMode())); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
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
             std::cout << "... m svfit : " << algo.mass() << " +/- " << algo.massUncert() << std::endl;
          } else {
            std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
          }

          delete inputFile_visPtResolution;

          return svfitResult;
    }

//    void FillSyncTree(const edm::Event& iEvent, const analysis::SelectionResultsV2& selection, Run2::SyncTree syncTree)
//        {
//            using namespace analysis;
//            static const float default_value = Run2::DefaultFloatFillValueForSyncTree();

//            // Event
//            std::cout<<"~~~~~~~~~~~~~EVENT Info~~~~~~~~~"<<std::endl;
//            std::cout<<"Run = "<<iEvent.id().run()<<"  Lumi = "<<iEvent.id().luminosityBlock()
//                    <<" Event = "<<iEvent.id().event()<<std::endl;
//            syncTree.run()  = iEvent.id().run();
//            syncTree.lumi() = iEvent.id().luminosityBlock();
//            syncTree.evt()  = iEvent.id().event();

//            if(computeHT_){
//                syncTree.HT()   = selection.HT;
//                if(selection.HT<100) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::lt100);
//                if(100<=selection.HT && selection.HT<200) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f100to200);
//                if(200<=selection.HT && selection.HT<400) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f200to400);
//                if(400<=selection.HT && selection.HT<600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f400to600);
//                if(selection.HT>=600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::gt600);
//    //        syncTree->eventType() = static_cast<int>(selection.eventType);
//    //        syncTree->eventEnergyScale() = static_cast<int>(eventEnergyScale);
//            }
//            else {
//                syncTree.HT() = default_value;
//                syncTree.HTBin() = -1;
//            }
//    }


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

    analysis::CandidateV2PtrVector CollectBJets()
    {
        return CollectCandidateObjects("bjets", &BaseEDAnalyzer::SelectBJets, *jets);
    }

    virtual void SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
            throw std::runtime_error("BJets selection not implemented");
    }




};

    // CandidatePtrVector CollectSignalMuons()
    // {
    //     return CollectCandidateObjects("muons_sgn", &BaseAnalyzer::SelectSignalMuon, event->muons());
    // }
    //
    // virtual void SelectMuon(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Muon selection for signal not implemented");
    // }
    //
    // virtual void SelectSignalMuon(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Signal muon selection for signal not implemented");
    // }
    //
    // CandidatePtrVector CollectBackgroundMuons()
    // {
    //     return CollectCandidateObjects("muons_bkg", &BaseAnalyzer::SelectBackgroundMuon, event->muons());
    // }
    //
    // virtual void SelectBackgroundMuon(const CandidatePtr& muon, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     using namespace cuts::Htautau_Summer13::muonVeto;
    //     const ntuple::Muon& object = muon->GetNtupleObject<ntuple::Muon>();
    //
    //     cut(true, ">0 mu cand");
    //     cut(X(pt) > pt, "pt");
    //     cut(std::abs( X(eta) ) < eta, "eta");
    //     const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
    //     cut(Y(DeltaZ)  < dz, "dz");
    //     const TVector3 mu_vertex(object.vx, object.vy, object.vz);
    //     const double d0_PV = Calculate_dxy(mu_vertex, primaryVertex->GetPosition(), muon->GetMomentum());
    //     cut(std::abs( Y(d0_PV) ) < d0, "d0");
    //     cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
    //     cut(X(isPFMuon) == isPFMuon, "PF");
    //     cut(X(nMatchedStations) > nMatched_Stations, "stations");
    //     cut(X(pixHits) > pixHits, "pix_hits");
    //     cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
    //     cut(X(pfRelIso) < pfRelIso, "pFRelIso");
    // }
    //
    // CandidatePtrVector CollectTaus()
    // {
    //     return CollectCandidateObjects("taus", &BaseAnalyzer::SelectTau, correctedTaus);
    // }
    //
    // CandidatePtrVector CollectSignalTaus()
    // {
    //     return CollectCandidateObjects("taus_sgn", &BaseAnalyzer::SelectSignalTau, correctedTaus);
    // }
    //
    // virtual void SelectTau(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Tau selection for signal not implemented");
    // }
    //
    // virtual void SelectSignalTau(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Signal tau selection for signal not implemented");
    // }
    //
    // CandidatePtrVector CollectElectrons()
    // {
    //     return CollectCandidateObjects("electrons", &BaseAnalyzer::SelectElectron, event->electrons());
    // }
    //
    // CandidatePtrVector CollectSignalElectrons()
    // {
    //     return CollectCandidateObjects("electrons_sgn", &BaseAnalyzer::SelectSignalElectron, event->electrons());
    // }
    //
    // virtual void SelectElectron(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Electron selection for signal not implemented");
    // }
    //
    // virtual void SelectSignalElectron(const CandidatePtr& candidate, SelectionManager& selectionManager,
    //                                   cuts::Cutter& cut)
    // {
    //     throw std::runtime_error("Electron selection for signal not implemented");
    // }
    //
    // CandidatePtrVector CollectBackgroundElectrons()
    // {
    //     return CollectCandidateObjects("electrons_bkg", &BaseAnalyzer::SelectBackgroundElectron, event->electrons());
    // }

//     virtual void SelectBackgroundElectron(const CandidatePtr& electron, SelectionManager& selectionManager,
//                                           cuts::Cutter& cut)
//     {
//         using namespace cuts::Htautau_Summer13::electronVeto;
//         const ntuple::Electron& object = electron->GetNtupleObject<ntuple::Electron>();
//
//         cut(true, ">0 ele cand");
//         cut(X(pt) > pt, "pt");
//         const double eta = std::abs( X(eta) );
//         cut(eta < eta_high, "eta");
//         const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
//         cut(Y(DeltaZ)  < dz, "dz");
//         const TVector3 ele_vertex(object.vx, object.vy, object.vz);
//         // same as dB
//         const double d0_PV = Calculate_dxy(ele_vertex, primaryVertex->GetPosition(), electron->GetMomentum());
//         cut(std::abs( Y(d0_PV) ) < d0, "d0");
//         cut(X(pfRelIso) < pFRelIso, "pFRelIso");
//         const size_t pt_index = object.pt < ref_pt ? 0 : 1;
//         const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
//         cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
//         cut(X(missingHits) < missingHits, "missingHits");
//         cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
//     }
//
//     CandidatePtrVector CollectBJets(const CandidatePtrVector& looseJets, bool doReTag, bool applyCsvCut)
//     {
//         using namespace cuts::Htautau_Summer13::btag;
//
//         static const std::map<EventEnergyScale, std::pair<int,int>> btag_modes_map = {
//             { EventEnergyScale::Central, {0,0} }, { EventEnergyScale::TauUp, {0,0} },
//             { EventEnergyScale::TauDown, {0,0} }, { EventEnergyScale::JetUp, {0,0} },
//             { EventEnergyScale::JetDown, {0,0} }, { EventEnergyScale::BtagEfficiencyUp, {2,0} },
//             { EventEnergyScale::BtagEfficiencyDown, {1,0} }, { EventEnergyScale::BtagFakeUp, {0,2} },
//             { EventEnergyScale::BtagFakeDown, {0,1} }
//         };
//
//         CandidatePtrVector bjets;
//         const std::pair<int,int> btag_pair = btag_modes_map.at(eventEnergyScale);
//         const int btag_mode = btag_pair.first;
//         const int bfake_mode = btag_pair.second;
//         for(const CandidatePtr& looseJetCandidate : looseJets) {
//             const ntuple::Jet& looseJet = looseJetCandidate->GetNtupleObject<ntuple::Jet>();
//             if(looseJet.pt <= pt || std::abs(looseJet.eta) >= eta) continue;
//             if(doReTag && !btag::ReTag(looseJet, btag::payload::EPS13, btag::tagger::CSVM, btag_mode, bfake_mode, CSV))
//                 continue;
//             else if(!doReTag && applyCsvCut && looseJet.combinedSecondaryVertexBJetTags <= CSV)
//                 continue;
//
//             bjets.push_back(looseJetCandidate);
//         }
//
//         const auto bjetsSelector = [&] (const CandidatePtr& first, const CandidatePtr& second) -> bool
//         {
//             const ntuple::Jet& first_bjet = first->GetNtupleObject<ntuple::Jet>();
//             const ntuple::Jet& second_bjet = second->GetNtupleObject<ntuple::Jet>();
//
//             return first_bjet.combinedSecondaryVertexBJetTags > second_bjet.combinedSecondaryVertexBJetTags;
//         };
//
//         std::sort(bjets.begin(), bjets.end(), bjetsSelector);
//         return bjets;
//     }
//
//     CandidatePtrVector CollectLooseJets()
//     {
//         return CollectCandidateObjects("loose_jets", &BaseAnalyzer::SelectLooseJet, GetNtupleJets());
//     }
//
//     virtual void SelectLooseJet(const CandidatePtr& jet, SelectionManager& selectionManager, cuts::Cutter& cut)
//     {
//         using namespace cuts::Htautau_Summer13::jetID;
//         const ntuple::Jet& object = jet->GetNtupleObject<ntuple::Jet>();
//
//         cut(true, ">0 jet cand");
//         cut(X(pt) > pt_loose, "pt");
//         cut(std::abs( X(eta) ) < eta, "eta");
//         const bool passLooseID = passPFLooseId(object);
//         cut(Y(passLooseID) == pfLooseID, "pfLooseID");
//         const bool passPUlooseID = ntuple::JetID_MVA::PassLooseId(object.puIdBits);
//         cut(Y(passPUlooseID) == puLooseID, "puLooseID");
//     }
//
//     CandidatePtrVector CollectJets(const CandidatePtrVector& looseJets)
//     {
//         using namespace cuts::Htautau_Summer13::jetID;
//         CandidatePtrVector jets;
//         for(const CandidatePtr& looseJet : looseJets) {
//             if(looseJet->GetMomentum().Pt() > pt)
//                 jets.push_back(looseJet);
//         }
//
//         const auto jetsSelector = [&] (const CandidatePtr& first_jet, const CandidatePtr& second_jet) -> bool
//         {
//             const double first_pt = first_jet->GetMomentum().Pt();
//             const double second_pt = second_jet->GetMomentum().Pt();
//             return first_pt > second_pt;
//         };
//
//         std::sort(jets.begin(), jets.end(), jetsSelector);
//
//         return jets;
//     }
//
//     VertexPtrVector CollectVertices()
//     {
//         const auto base_selector = [&](const VertexPtr& vertex, SelectionManager& selectionManager,
//                                        cuts::Cutter& cut)
//             { SelectVertex(vertex, selectionManager, cut); };
//         const auto vertex_comparitor = [&](const VertexPtr& first, const VertexPtr& second) -> bool
//             { return *first < *second; };
//         return CollectObjects<Vertex>("vertices", base_selector, event->vertices(), vertex_comparitor);
//     }
//
//     void SelectVertex(const VertexPtr& vertex, SelectionManager& selectionManager, cuts::Cutter& cut)
//     {
//         using namespace cuts::Htautau_Summer13::vertex;
//         const ntuple::Vertex& object = vertex->GetNtupleObject();
//
//         cut(true, ">0 vertex");
//         cut(X(ndf) > ndf, "ndf");
//         cut(std::abs( X(z) ) < z, "z");
//         const double r_vertex = std::sqrt(object.x*object.x+object.y*object.y);
//         cut(std::abs( Y(r_vertex) ) < r, "r");
//     }
//
//     CandidatePtrVector FindCompatibleObjects(const CandidatePtrVector& objects1, const CandidatePtrVector& objects2,
//                                           double minDeltaR, Candidate::Type type, const std::string& hist_name,
//                                           int expectedCharge = Candidate::UnknownCharge())
//     {
//         CandidatePtrVector result;
//         for(const CandidatePtr& object1 : objects1) {
//             for(const CandidatePtr& object2 : objects2) {
//                 if(object2->GetMomentum().DeltaR(object1->GetMomentum()) > minDeltaR) {
//                     const CandidatePtr candidate(new Candidate(type, object1, object2));
//                     if (expectedCharge != Candidate::UnknownCharge() && candidate->GetCharge() != expectedCharge)
//                         continue;
//                     result.push_back(candidate);
//                     GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),
//                                                       GetEventWeights().GetPartialWeight());
//                 }
//             }
//         }
//         GetAnaData().N_objects(hist_name).Fill(result.size(), GetEventWeights().GetPartialWeight());
//         return result;
//     }
//
//
//     CandidatePtrVector FindCompatibleObjects(const CandidatePtrVector& objects, double minDeltaR, Candidate::Type type,
//                                           const std::string& hist_name, int expectedCharge = Candidate::UnknownCharge())
//     {
//         CandidatePtrVector result;
//         for (unsigned n = 0; n < objects.size(); ++n){
//             for (unsigned k = n+1; k < objects.size(); ++k){
// //                std::cout << "first tau momentum " << objects.at(n).momentum << std::endl;
// //                std::cout << "second tau momentum " << objects.at(k).momentum << std::endl;
// //                std::cout << "DeltaR " << objects.at(n).momentum.DeltaR(objects.at(k).momentum) << std::endl;
// //                std::cout << "first tau charge " << objects.at(n).charge << std::endl;
// //                std::cout << "second tau charge " << objects.at(k).charge << std::endl;
//                 if(objects.at(n)->GetMomentum().DeltaR(objects.at(k)->GetMomentum()) > minDeltaR) {
//                     const CandidatePtr candidate(new Candidate(type, objects.at(n), objects.at(k)));
//                     if (expectedCharge != Candidate::UnknownCharge() && candidate->GetCharge() != expectedCharge )
//                         continue;
//                     result.push_back(candidate);
//                     GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),
//                                                       GetEventWeights().GetPartialWeight());
//                 }
//             }
//         }
//         GetAnaData().N_objects(hist_name).Fill(result.size(), GetEventWeights().GetPartialWeight());
//         return result;
//     }
//
//     CandidatePtrVector FilterCompatibleObjects(const CandidatePtrVector& objectsToFilter,
//                                                const CandidatePtr& referenceObject, double minDeltaR)
//     {
//         CandidatePtrVector result;
//         for(const CandidatePtr& filterObject : objectsToFilter) {
//             bool allDaughterPassed = true;
//             for (const CandidatePtr& daughter : referenceObject->GetFinalStateDaughters()){
//                 if(filterObject->GetMomentum().DeltaR(daughter->GetMomentum()) <= minDeltaR) {
//                     allDaughterPassed = false;
//                     break;
//                 }
//             }
//             if (allDaughterPassed) result.push_back(filterObject);
//         }
//         return result;
//     }
//
//     ntuple::TauVector ApplyTauCorrections(const VisibleGenObjectVector& hadronic_taus, bool useLegacyCorrections)
//     {
//         using namespace cuts::Htautau_Summer13::tauCorrections;
//
//         ntuple::TauVector correctedTaus;
//
//         for(const ntuple::Tau& tau : GetNtupleTaus()) {
//             TLorentzVector momentum;
//             momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
//
//             const bool hasMCmatch = FindMatchedObjects(momentum, hadronic_taus, deltaR).size() != 0;
//             const double scaleFactor = MomentumScaleFactor(hasMCmatch, momentum.Pt(),
//                                    ntuple::tau_id::ConvertToHadronicDecayMode(tau.decayMode), useLegacyCorrections);
//             const TLorentzVector correctedMomentum = momentum * scaleFactor;
//             ntuple::Tau correctedTau(tau);
//             correctedTau.pt = correctedMomentum.Pt();
//             correctedTau.eta = correctedMomentum.Eta();
//             correctedTau.phi = correctedMomentum.Phi();
//             correctedTau.mass = correctedMomentum.M();
//             correctedTaus.push_back(correctedTau);
//         }
//         return correctedTaus;
//     }
//
//     ntuple::MET ApplyTauCorrectionsToMVAMET(const ntuple::MET& metMVA, const ntuple::TauVector& correctedTaus)
//     {
//         TLorentzVector sumCorrectedTaus, sumTaus;
//         for(const ntuple::Tau& tau : GetNtupleTaus()) {
//             TLorentzVector momentum;
//             momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
//             sumTaus += momentum;
//         }
//
//         for (const ntuple::Tau& correctedTau : correctedTaus) {
//             TLorentzVector correctedMomentum;
//             correctedMomentum.SetPtEtaPhiM(correctedTau.pt, correctedTau.eta, correctedTau.phi,correctedTau.mass);
//             sumCorrectedTaus += correctedMomentum;
//         }
//
//         TLorentzVector met, metCorrected;
//         met.SetPtEtaPhiM(metMVA.pt, 0, metMVA.phi, 0.);
//         metCorrected = met + sumTaus - sumCorrectedTaus;
//         ntuple::MET correctedMET = metMVA;
//         correctedMET.pt = metCorrected.Pt();
//         correctedMET.phi = metCorrected.Phi();
//         return correctedMET;
//     }
//
//     CandidatePtrVector ApplyTriggerMatch(const CandidatePtrVector& higgses, const std::set<std::string>& hltPaths,
//                                          bool useStandardTriggerMatch)
//     {
//         CandidatePtrVector triggeredHiggses;
//         for (const auto& higgs : higgses){
//             if(!useStandardTriggerMatch && HaveTriggerMatched(event->triggerObjects(), hltPaths, *higgs,
//                                                                         cuts::Htautau_Summer13::DeltaR_triggerMatch))
//                 triggeredHiggses.push_back(higgs);
//             if (useStandardTriggerMatch && HaveTriggerMatched(hltPaths, *higgs))
//                 triggeredHiggses.push_back(higgs);
//         }
//         return triggeredHiggses;
//     }
//
//     ntuple::MET ApplyRecoilCorrections(const CandidatePtr& higgs, const GenParticle* resonance,
//                                        const size_t njets, const ntuple::MET& correctedMET)
//     {
//         if (config.ApplyRecoilCorrection()){
//             if(!resonance && config.ApplyRecoilCorrectionForW())
//                 resonance = FindWboson();
//             if(!resonance && config.ApplyRecoilCorrectionForZ()){
//                 GenParticlePtrVector ZProducts;
//                 bool ztt;
//                 resonance = FindZboson(ZProducts,ztt);
//             }
//             if(resonance)
//                 return GetRecoilCorrectionProducer().ApplyCorrection(correctedMET, higgs->GetMomentum(),
//                                                                      resonance->momentum, njets);
//         }
//         return correctedMET;
//     }
//
//     const GenParticle* FindWboson()
//     {
//         static const particles::ParticleCodes Wcode = { particles::W_plus };
//         static const particles::ParticleCodes WDecay_tau = { particles::tau, particles::nu_tau };
//         static const particles::ParticleCodes WDecay_electron = { particles::e, particles::nu_e };
//         static const particles::ParticleCodes WDecay_muon = { particles::mu, particles::nu_mu };
//
//         const GenParticleSet Wparticles_all = genEvent.GetParticles(Wcode);
//
//         GenParticleSet Wparticles;
//         for(const GenParticle* w : Wparticles_all) {
//             if(w->mothers.size() == 1) {
//                 const GenParticle* mother = w->mothers.at(0);
//                 if(mother->pdg.Code == particles::W_plus && mother->status == particles::HardInteractionProduct)
//                     Wparticles.insert(w);
//             }
//          }
//
//         if (Wparticles.size() == 0) return nullptr;
//
//         if (Wparticles.size() > 1)
//             throw exception("more than 1 W in the event");
//
//         const GenParticle* Wboson = *Wparticles.begin();
//         while(Wboson->daughters.size() == 1 && Wboson->daughters.front()->pdg.Code == particles::W_plus)
//             Wboson = Wboson->daughters.front();
//
//         GenParticlePtrVector WProducts;
//         if(FindDecayProducts(*Wboson, WDecay_tau, WProducts,true) ||
//                 FindDecayProducts(*Wboson, WDecay_electron, WProducts,true) ||
//                 FindDecayProducts(*Wboson, WDecay_muon, WProducts,true))
//             return Wboson;
//
//         throw exception("not leptonic W decay");
//
//     }
//
//     const GenParticle* FindZboson(GenParticlePtrVector& ZProducts, bool& ztt)
//     {
//         static const particles::ParticleCodes Zcode = { particles::Z };
//         static const particles::ParticleCodes ZDecay_electrons = { particles::e, particles::e };
//         static const particles::ParticleCodes ZDecay_muons = { particles::mu, particles::mu };
//         static const particles::ParticleCodes ZDecay_taus = { particles::tau, particles::tau };
//
//         const GenParticleSet Zparticles_all = genEvent.GetParticles(Zcode);
//
//         GenParticleSet Zparticles;
//         for(const GenParticle* z : Zparticles_all) {
//             const bool is_hard_interaction_z = z->mothers.size() == 1 && z->mothers.front()->pdg.Code == particles::Z
//                     && z->mothers.front()->status == particles::HardInteractionProduct;
//             const bool is_pp_z = z->mothers.size() == 2 && z->mothers.at(0)->pdg.Code == particles::p
//                     && z->mothers.at(1)->pdg.Code == particles::p;
//             if(is_hard_interaction_z || is_pp_z)
//                 Zparticles.insert(z);
//          }
//
//         if (Zparticles.size() > 1 || Zparticles.size() == 0)
//             throw exception("not 1 Z per event");
//
//         const GenParticle* Z_mc = *Zparticles.begin();
//         while(Z_mc->daughters.size() == 1 && Z_mc->daughters.front()->pdg.Code == particles::Z)
//             Z_mc = Z_mc->daughters.front();
//
//         ztt = FindDecayProducts(*Z_mc, ZDecay_taus, ZProducts, true);
//         if (ztt || FindDecayProducts(*Z_mc, ZDecay_electrons, ZProducts, true)
//                  || FindDecayProducts(*Z_mc, ZDecay_muons, ZProducts, true))
//             return Z_mc;
//         throw exception("not leptonic Z decay");
//     }
//
//     CandidatePtr SelectSemiLeptonicHiggs(const CandidatePtrVector& higgses)
//     {
//         if(!higgses.size())
//             throw std::runtime_error("no available higgs candidate to select");
//         const auto higgsSelector = [&] (const CandidatePtr& first, const CandidatePtr& second) -> bool
//         {
//             const double first_Pt1 = first->GetDaughters().at(0)->GetMomentum().Pt();
//             const double first_Pt2 = first->GetDaughters().at(1)->GetMomentum().Pt();
//             const double first_sumPt = first_Pt1 + first_Pt2;
//             const double second_Pt1 = second->GetDaughters().at(0)->GetMomentum().Pt();
//             const double second_Pt2 = second->GetDaughters().at(1)->GetMomentum().Pt();
//             const double second_sumPt = second_Pt1 + second_Pt2;
//
//             return first_sumPt < second_sumPt;
//         };
//         return *std::max_element(higgses.begin(), higgses.end(), higgsSelector) ;
//     }
//
//     bool FindAnalysisFinalState(finalState::bbTauTau& final_state)
//     {
//         static const particles::ParticleCodes resonanceCodes = { particles::MSSM_H, particles::MSSM_A };
//         static const particles::ParticleCodes resonanceDecay_1 = { particles::Higgs, particles::Higgs };
//         static const particles::ParticleCodes resonanceDecay_2 = { particles::Z, particles::Higgs };
//         static const particles::ParticleCodes SM_ResonanceCodes = { particles::Higgs, particles::Z,
//                                                                     particles::MSSM_H, particles::MSSM_A };
//         static const particles::ParticleCodes SM_ResonanceDecay_1 = { particles::tau, particles::tau };
//         static const particles::ParticleCodes SM_ResonanceDecay_2 = { particles::b, particles::b };
//         static const particles::ParticleCodes2D HiggsDecays = { SM_ResonanceDecay_1, SM_ResonanceDecay_2 };
//
//         genEvent.Initialize(event->genParticles());
//         final_state.Reset();
//
//         const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
//
//         if (resonances.size() > 1)
//             throw exception("more than 1 resonance per event");
//
//         if (resonances.size() == 1) {
//             final_state.resonance = *resonances.begin();
//
//             bool doubleHiggsSignal = true;
//             GenParticlePtrVector HiggsBosons;
//             if(!FindDecayProducts(*final_state.resonance, resonanceDecay_1, HiggsBosons) &&
//                     !FindDecayProducts(*final_state.resonance, resonanceDecay_2, HiggsBosons))
//                 doubleHiggsSignal = false;
//
//             if(doubleHiggsSignal) {
//                 GenParticleVector2D HiggsDecayProducts;
//                 GenParticleIndexVector HiggsIndexes;
//                 GenParticlePtrVector Higgs_ToTauTau_Product;
//                 GenParticlePtrVector Higgs_ToBB_Product;
//                 const bool HH_bbtautau = FindDecayProducts2D(HiggsBosons, HiggsDecays, HiggsDecayProducts, HiggsIndexes);
//                 if(HH_bbtautau){
//                     Higgs_ToTauTau_Product = HiggsDecayProducts.at(0);
//                     Higgs_ToBB_Product = HiggsDecayProducts.at(1);
//                     final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(0));
//                     final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(1));
//                 } else if (FindDecayProducts(*HiggsBosons.at(1),SM_ResonanceDecay_1,Higgs_ToTauTau_Product)){
//                     final_state.Higgs_TauTau = HiggsBosons.at(1);
//                 } else if (FindDecayProducts(*HiggsBosons.at(1),SM_ResonanceDecay_2,Higgs_ToBB_Product)){
//                     final_state.Higgs_BB = HiggsBosons.at(1);
//                 } else
//                     throw exception("Nor HH-> bbtautau, nor A->Zh->lltautau, nor A->Zh->llbb");
//
//                 for(const GenParticle* tau : Higgs_ToTauTau_Product) {
//                     const VisibleGenObject tau_products(tau);
//                     final_state.taus.push_back(tau_products);
// //                    if(tau_products.finalStateChargedHadrons.size() != 0)
//                     if(!IsLeptonicTau(*tau)){
//
//                         final_state.hadronic_taus.push_back(tau_products);
//                     }
//                 }
//                 for(const GenParticle* b : Higgs_ToBB_Product)
//                     final_state.b_jets.push_back(VisibleGenObject(b));
//
//                 return HH_bbtautau;
//             }
//         }
//
//         if(config.ExpectedOneNonSMResonance())
//             throw exception("Non-SM resonance not found.");
//
//         //search H->bb, H->tautau
//         const GenParticleSet SM_particles = genEvent.GetParticles(SM_ResonanceCodes);
//
//         GenParticlePtrVector SM_ResonanceToTauTau_products;
//         GenParticlePtrVector SM_ResonanceToBB_products;
//
//         for (const GenParticle* SM_particle : SM_particles){
//             GenParticlePtrVector resonanceDecayProducts;
//             if(FindDecayProducts(*SM_particle, SM_ResonanceDecay_1,resonanceDecayProducts)){
//                 if(!final_state.Higgs_TauTau || (final_state.Higgs_TauTau->pdg.Code != particles::Higgs
//                                              && SM_particle->pdg.Code == particles::Higgs)) {
//                     final_state.Higgs_TauTau = SM_particle;
//                     SM_ResonanceToTauTau_products = resonanceDecayProducts;
//                 } else if((final_state.Higgs_TauTau->pdg.Code == particles::Higgs
//                            && SM_particle->pdg.Code == particles::Higgs)
//                           || (final_state.Higgs_TauTau->pdg.Code != particles::Higgs
//                               && SM_particle->pdg.Code != particles::Higgs)) {
//                     throw exception("more than one SM resonance to tautau per event");
//                 }
//             }
//             else if (FindDecayProducts(*SM_particle, SM_ResonanceDecay_2,resonanceDecayProducts)){
//                 if(!final_state.Higgs_BB || (final_state.Higgs_BB->pdg.Code != particles::Higgs
//                                              && SM_particle->pdg.Code == particles::Higgs)) {
//                     final_state.Higgs_BB = SM_particle;
//                     SM_ResonanceToBB_products = resonanceDecayProducts;
//                 } else if((final_state.Higgs_BB->pdg.Code == particles::Higgs
//                            && SM_particle->pdg.Code == particles::Higgs)
//                           || (final_state.Higgs_BB->pdg.Code != particles::Higgs
//                               && SM_particle->pdg.Code != particles::Higgs)) {
//                     throw exception("more than one SM resonance to bb per event");
//                 }
//             }
//         }
//
//         for(const GenParticle* tau : SM_ResonanceToTauTau_products) {
//             const VisibleGenObject tau_products(tau);
//             final_state.taus.push_back(tau_products);
//             if(!IsLeptonicTau(*tau))
//                 final_state.hadronic_taus.push_back(tau_products);
//         }
//
//         for(const GenParticle* b : SM_ResonanceToBB_products)
//             final_state.b_jets.push_back(VisibleGenObject(b));
//
//         if (!final_state.Higgs_TauTau && !final_state.Higgs_BB) {
//             if(config.ExpectedAtLeastOneSMResonanceToTauTauOrToBB())
//                 throw exception("SM resonance to tautau or to bb not found.");
//             return false;
//         }
//
//         return true;
//     }
//
//     ntuple::EventType DoEventCategorization(const CandidatePtr& higgs)
//     {
//         using namespace cuts::Htautau_Summer13::DrellYannCategorization;
//         if(!config.DoZEventCategorization())
//             return ntuple::EventType::Unknown;
//
//         GenParticlePtrVector ZProducts;
//         bool ztt;
//         FindZboson(ZProducts,ztt);
//
//         static const particles::ParticleCodes light_lepton_codes = { particles::e, particles::mu };
//
//         const GenParticleSet light_leptons = genEvent.GetParticles(light_lepton_codes, minimal_genParticle_pt);
//         const CandidatePtrVector hadronic_taus = higgs->GetDaughters(Candidate::Type::Tau);
//
//         size_t n_hadronic_matches = 0, n_leptonic_matches = 0;
//         for(const CandidatePtr& reco_tau : hadronic_taus) {
//
//             for(const GenParticle* gen_product : ZProducts) {
//                 const VisibleGenObject visible_gen_object(gen_product);
// //                 std::cout <<  "GenVisibleTau: " << visible_gen_object.visibleMomentum <<
// //                              "; NofLeptons: " << visible_gen_object.finalStateChargedLeptons.size() <<
// //                             "; GenTauOrigin: " << visible_gen_object.origin->momentum <<
// //                               "; pdg: " << visible_gen_object.origin->pdg.Code.Name() << ", status= " <<
// //                               visible_gen_object.origin->status << std::endl;
//                 if(gen_product->pdg.Code != particles::tau || IsLeptonicTau(*gen_product) ||
//                         visible_gen_object.visibleMomentum.Pt() <= minimal_visible_momentum) continue;
//                 if(HasMatchWithMCObject(reco_tau->GetMomentum(), &visible_gen_object, deltaR_matchGenParticle, true)) {
//                     ++n_hadronic_matches;
//                     break;
//                 }
//             }
//
//             for(const GenParticle* gen_product : light_leptons) {
//                 if(HasMatchWithMCParticle(reco_tau->GetMomentum(), gen_product, deltaR_matchGenParticle)) {
//                     ++n_leptonic_matches;
//                     break;
//                 }
//             }
//         }
//         //genEvent.Print();
//
//         if(ztt && n_hadronic_matches == hadronic_taus.size()) return ntuple::EventType::ZTT;
//         if(n_leptonic_matches) return ztt ? ntuple::EventType::ZTT_L : ntuple::EventType::ZL;
//         return ntuple::EventType::ZJ;
//     }
//
//     bool GenFilterForZevents(const finalState::bbTauTau& final_state)
//     {
//         using namespace cuts::Htautau_Summer13::DYEmbedded;
//         if (final_state.taus.size() != 2)
//             throw exception("not 2 taus in the event at Gen Level");
//         const GenParticle* firstTau = final_state.taus.at(0).origin;
//         const GenParticle* secondTau = final_state.taus.at(1).origin;
//         if ((firstTau->momentum + secondTau->momentum).M() > invariantMassCut) return true;
//         return false;
//     }
//
//     kinematic_fit::four_body::FitResults RunKinematicFit(const CandidatePtrVector& bjets,
//                                                          const Candidate& higgs_to_taus, const ntuple::MET& met)
//     {
//         using namespace kinematic_fit::four_body;
//         using namespace cuts::Htautau_Summer13;
//
//         if(bjets.size() < 2)
//             return FitResults();
//
//         const TLorentzVector met_momentum = MakeLorentzVectorPtEtaPhiM(met.pt, 0, met.phi, 0);
//         const TMatrix met_cov = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);
//
//         const FitInput input(bjets.at(0)->GetMomentum(), bjets.at(1)->GetMomentum(),
//                              higgs_to_taus.GetDaughters().at(0)->GetMomentum(),
//                              higgs_to_taus.GetDaughters().at(1)->GetMomentum(),
//                              met_momentum, met_cov);
//         return Fit(input);
//     }
//
//     ntuple::MET ComputeMvaMet(const CandidatePtr& higgs, const VertexPtrVector& goodVertices)
//     {
//         CandidatePtrVector originalDaughters;
//         for(const CandidatePtr& daughter : higgs->GetDaughters()) {
//             CandidatePtr original = daughter;
//             if(daughter->GetType() == Candidate::Type::Tau) {
//                 const ntuple::Tau* ntuple_tau = &daughter->GetNtupleObject<ntuple::Tau>();
//                 const auto position = ntuple_tau - &(*correctedTaus.begin());
//                 original = CandidatePtr(new Candidate(GetNtupleTaus().at(position)));
//             }
//             originalDaughters.push_back(original);
//         }
//         CandidatePtr originalHiggs(new Candidate(higgs->GetType(), originalDaughters.at(0), originalDaughters.at(1)));
//         return mvaMetProducer.ComputeMvaMet(originalHiggs, event->pfCandidates(), GetNtupleJets(), primaryVertex,
//                                             goodVertices);
//     }

// protected:
//     Config config;
//     std::shared_ptr<tools::ProgressReporter> progressReporter;
//     std::shared_ptr<const EventDescriptor> event;
//     std::shared_ptr<TreeExtractor> treeExtractor;
//     std::shared_ptr<TFile> outputFile;
//     root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
//     size_t maxNumberOfEvents;
//     GenEvent genEvent;
//     VertexPtr primaryVertex;
//     MvaMetProducer mvaMetProducer;
//     ntuple::TauVector correctedTaus;
//     EventEnergyScale eventEnergyScale;
//     std::shared_ptr<JetEnergyUncertaintyCorrector> jetEnergyUncertaintyCorrector;
//
// private:
//     ntuple::TauVector scaledTaus;
//     ntuple::JetVector scaledJets;
