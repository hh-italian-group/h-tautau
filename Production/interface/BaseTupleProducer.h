/*! Definition of BaseTupleProducer class: the base class for all X->HH->bbTauTau and H->tautau tuple producers.
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//Trigger

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
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

//HHbbTauTau Framework
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/CutTools.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/Candidate.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "h-tautau/Cuts/include/H_tautau_2016_mssm.h"
#include "h-tautau/Cuts/include/H_tautau_2016_sm.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/EventLoader.h"

//SVFit
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "SelectionResults.h"

//Recoil Correction
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "TriggerTools.h"

struct TupleProducerData : root_ext::AnalyzerData {
    using AnalyzerData::AnalyzerData;
    SELECTION_ENTRY(Selection)
    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};

namespace analysis {
namespace detail {
template<>
inline bool CompareIsolations<pat::Electron>(double iso_1, double iso_2) { return iso_1 < iso_2; }
template<>
inline bool CompareIsolations<pat::Muon>(double iso_1, double iso_2) { return iso_1 < iso_2; }
template<>
inline bool CompareIsolations<pat::Tau>(double iso_1, double iso_2) { return iso_1 > iso_2; }
}
}

enum class ProductionMode { hh, h_tt, h_tt_mssm,h_tt_sm };

ENUM_NAMES(ProductionMode) = {
    { ProductionMode::hh, "hh" },
    { ProductionMode::h_tt, "h_tt" },
    { ProductionMode::h_tt_mssm, "h_tt_mssm" },
    { ProductionMode::h_tt_sm, "h_tt_sm" },
};

class BaseTupleProducer : public edm::EDAnalyzer {
public:
    using ElectronCandidate = analysis::LeptonCandidate<pat::Electron, edm::Ptr<pat::Electron>>;
    using MuonCandidate = analysis::LeptonCandidate<pat::Muon>;
    using TauCandidate = analysis::LeptonCandidate<pat::Tau>;
    using JetCandidate = analysis::Candidate<pat::Jet>;
    using MET = analysis::MissingET<pat::MET>;
    using MetCovMatrix = MET::CovMatrix;
    using SelectionManager = analysis::SelectionManager;
    using Cutter = cuts::Cutter<SelectionManager>;
    using LorentzVector = analysis::LorentzVector;
    using LorentzVectorM = analysis::LorentzVectorM;
    using LorentzVectorE = analysis::LorentzVectorE;

private:
    TupleProducerData anaData;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
    edm::EDGetToken electronsMiniAOD_token;
    edm::EDGetTokenT<edm::ValueMap<bool>> eleTightIdMap_token, eleMediumIdMap_token, eleCutBasedVetoMap_token;
    edm::EDGetToken tausMiniAOD_token;
    edm::EDGetToken muonsMiniAOD_token;
    edm::EDGetToken vtxMiniAOD_token;
    edm::EDGetToken pfMETAOD_token;
    edm::EDGetToken jetsMiniAOD_token;
    edm::EDGetToken fatJetsMiniAOD_token;
    edm::EDGetTokenT<MetCovMatrix> metCovMatrix_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUInfo_token;
    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_token;
    edm::EDGetTokenT<GenEventInfoProduct> genWeights_token;
    edm::EDGetTokenT<TtGenEvent> topGenEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<edm::View<reco::GenJet>> genJets_token;
    edm::EDGetTokenT<bool> badPFMuonFilter_token, badChCandidateFilter_token;

protected:
    const ProductionMode productionMode;
    const bool isMC, applyTriggerMatch, runSVfit, runKinFit, applyRecoilCorr;
    const int nJetsRecoilCorr;    
    const bool saveGenTopInfo, saveGenBosonInfo, saveGenJetInfo;
    std::vector<std::string> hltPaths;
    ntuple::EventTuple eventTuple;
    analysis::TriggerTools triggerTools;
    std::shared_ptr<analysis::sv_fit::FitProducer> svfitProducer;
    std::shared_ptr<analysis::kin_fit::FitProducer> kinfitProducer;
    std::shared_ptr<RecoilCorrector> recoilPFMetCorrector;

private:
    const edm::Event *edmEvent;
    edm::Handle<std::vector<pat::Electron> > pat_electrons;
    edm::Handle<std::vector<pat::Tau> > pat_taus;
    edm::Handle<std::vector<pat::Muon> > pat_muons;
    edm::Handle<edm::View<reco::Vertex> > vertices;
    edm::Handle<edm::View<pat::MET> > pfMETs;
    edm::Handle<std::vector<pat::Jet> > pat_jets;
    edm::Handle<std::vector<pat::Jet> > pat_fatJets;
    edm::Handle<MetCovMatrix> metCovMatrix;
    edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
    edm::Handle<LHEEventProduct> lheEventProduct;
    edm::Handle<GenEventInfoProduct> genEvt;
    edm::Handle<TtGenEvent> topGenEvent;

    edm::ESHandle<JetCorrectorParametersCollection> jetCorParColl;
    std::shared_ptr<JetCorrectionUncertainty> jecUnc;

    std::vector<analysis::EventEnergyScale> eventEnergyScales;

protected:
    edm::EventID eventId;
    analysis::EventEnergyScale eventEnergyScale;
    std::vector<ElectronCandidate> electrons;
    std::vector<MuonCandidate> muons;
    std::vector<TauCandidate> taus;
    std::vector<JetCandidate> jets, fatJets;
    std::shared_ptr<MET> met;
    edm::Ptr<reco::Vertex> primaryVertex;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    edm::Handle<edm::View<reco::GenJet>> genJets;
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions, medium_id_decisions, ele_cutBased_veto;

private:
    void InitializeAODCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void InitializeCandidateCollections(analysis::EventEnergyScale eventEnergyScale);
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    virtual void endJob() override;
    virtual void ProcessEvent(Cutter& cut) = 0;
    static double Isolation(const pat::Electron& electron);
    static double Isolation(const pat::Muon& muon);
    static double Isolation(const pat::Tau& tau);

public:
    BaseTupleProducer(const edm::ParameterSet& iConfig, const std::string& treeName);

protected:
    TupleProducerData& GetAnaData() { return anaData; }

    static bool PassPFLooseId(const pat::Jet& pat_jet);
    static bool PassICHEPMuonMediumId(const pat::Muon& pat_muon);

    double GetNumberOfPileUpInteractions() const;
    void ApplyBaseSelection(analysis::SelectionResultsBase& selection,
                            const std::vector<LorentzVector>& signalLeptonMomentums);
    void FillEventTuple(const analysis::SelectionResultsBase& selection,
                        const analysis::SelectionResultsBase* reference = nullptr);
    void FillElectronLeg(size_t leg_id, const ElectronCandidate& electron);
    void FillMuonLeg(size_t leg_id, const MuonCandidate& muon);
    void FillTauLeg(size_t leg_id, const TauCandidate& tau, bool fill_tauIds);
    void FillLheInfo(bool haveReference);
    void FillGenParticleInfo();
    void FillGenJetInfo();
    void FillLegGenMatch(size_t leg_id, const analysis::LorentzVectorXYZ& p4);
    void FillTauIds(size_t leg_id, const std::vector<pat::Tau::IdPair>& tauIds);
    void FillMetFilters();
    void ApplyRecoilCorrection(const std::vector<JetCandidate>& jets);

    std::vector<ElectronCandidate> CollectVetoElectrons(
            const std::vector<const ElectronCandidate*>& signalElectrons = {});
    std::vector<MuonCandidate> CollectVetoMuons(const std::vector<const MuonCandidate*>& signalMuons = {});
    std::vector<JetCandidate> CollectJets(const std::vector<LorentzVector>& signalLeptonMomentums);

    void SelectVetoElectron(const ElectronCandidate& electron, Cutter& cut,
                            const std::vector<const ElectronCandidate*>& signalElectrons) const;
    void SelectVetoMuon(const MuonCandidate& muon, Cutter& cut,
                        const std::vector<const MuonCandidate*>& signalMuons) const;
    void SelectJet(const JetCandidate& jet, Cutter& cut,
                   const std::vector<LorentzVector>& signalLeptonMomentums) const;

    template<typename Candidate1, typename Candidate2,
             typename ResultCandidate = analysis::CompositCandidate<Candidate1, Candidate2>>
    std::vector<ResultCandidate> FindCompatibleObjects(const std::vector<Candidate1>& objects1,
                                                       const std::vector<Candidate2>& objects2,
                                                       double minDeltaR, const std::string& hist_name,
                                                       int expectedCharge = analysis::AnalysisObject::UnknownCharge)
    {
        const double minDeltaR2 = std::pow(minDeltaR, 2);
        std::vector<ResultCandidate> result;
        for(const auto& object1 : objects1) {
            for(const auto& object2 : objects2) {
                if(ROOT::Math::VectorUtil::DeltaR2(object1.GetMomentum(), object2.GetMomentum()) > minDeltaR2) {
                    const ResultCandidate candidate(object1, object2);
                    if (expectedCharge !=analysis::AnalysisObject::UnknownCharge
                            && candidate.GetCharge() != expectedCharge) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.GetMomentum().M(), 1);
                }
            }
        }

        GetAnaData().N_objects(hist_name).Fill(result.size(), 1);
        return result;
    }

    template<typename Candidate, typename BaseSelectorType, typename Comparitor = std::less<Candidate>>
    std::vector<Candidate> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                              const std::vector<Candidate>& all_candidates,
                                              Comparitor comparitor = Comparitor())
    {
        static constexpr double weight = 1;
        std::ostringstream ss_suffix;
        ss_suffix << selection_label << "_" << eventEnergyScale;
        const std::string suffix = ss_suffix.str();
        cuts::ObjectSelector& objectSelector = GetAnaData().Selection(suffix);
        SelectionManager selectionManager(anaDataBeforeCut, suffix, weight);

        const auto selector = [&](size_t id) -> Candidate {
            const Candidate& candidate = all_candidates.at(id);

            Cutter cut(&objectSelector, &selectionManager);
            base_selector(candidate, cut);
            return candidate;
        };

        const auto selected = objectSelector.collect_objects<Candidate>(1, all_candidates.size(), selector, comparitor);

        SelectionManager selectionManager_afterCut(anaDataAfterCut, suffix, weight);
        for(const auto& candidate : selected) {
            Cutter cut(nullptr, &selectionManager_afterCut);
            base_selector(candidate, cut);
        }
        GetAnaData().N_objects(suffix).Fill(selected.size(), 1);
        GetAnaData().N_objects(suffix + "_original").Fill(all_candidates.size(), weight);

        return selected;
    }

    template<typename HiggsCandidate>
    static bool HiggsComparitor(const HiggsCandidate& h1, const HiggsCandidate& h2)
    {
        const auto& h1_leg1 = h1.GetFirstDaughter();
        const auto& h2_leg1 = h2.GetFirstDaughter();
        if(h1_leg1 != h2_leg1) {
            if(h1_leg1.GetIsolation() != h2_leg1.GetIsolation()) return h1_leg1.IsMoreIsolated(h2_leg1);
            if(h1_leg1->pt() != h2_leg1->pt()) return h1_leg1->pt() > h2_leg1->pt();
        }

        const auto& h1_leg2 = h1.GetSecondDaughter();
        const auto& h2_leg2 = h2.GetSecondDaughter();
        if(h1_leg2 != h2_leg2) {
            if(h1_leg2.GetIsolation() != h2_leg2.GetIsolation()) return h1_leg2.IsMoreIsolated(h2_leg2);
            if(h1_leg2->pt() != h2_leg2->pt()) return h1_leg2->pt() > h2_leg2->pt();
        }

        if(h1_leg1 == h2_leg1 && h1_leg2 == h2_leg2) return false;
        throw analysis::exception("not found a good criteria for best tau pair");
    }

    template<typename Candidate>
    static float GetUserFloat(const Candidate& obj, const std::string& name)
    {
        return obj->hasUserFloat(name) ? obj->userFloat(name) : ntuple::DefaultFillValue<float>();
    }
};
