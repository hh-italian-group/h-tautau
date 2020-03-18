/*! Definition of BaseTupleProducer class: the base class for all HH->bbTauTau and H->tautau tuple producers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Cuts/include/eleID_Run2.h"

#include "SelectionResults.h"
#include "TriggerTools.h"

struct TupleProducerData : public root_ext::AnalyzerData {
    explicit TupleProducerData(TDirectory* _directory, const std::string& subDirectoryName = "") :
        AnalyzerData(_directory, subDirectoryName)
    {
        Selection.SetMasterHist();
    }
    SELECTION_ENTRY(Selection)
    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};

struct SelectionData : public root_ext::AnalyzerData {
    explicit SelectionData(TDirectory* _directory, const std::string& subDirectoryName = "") :
        AnalyzerData(_directory, subDirectoryName) {}
    ANA_DATA_ENTRY(TH1D, h)
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

template<typename LeptonCandidate>
bool LeptonComparitor(const LeptonCandidate& l1, const LeptonCandidate& l2)
{
    return l1.IsMoreIsolated(l2);
}

template<>
inline bool LeptonComparitor(const SelectionResultsBase::ElectronCandidate& l1,
                             const SelectionResultsBase::ElectronCandidate& l2)
{
    static const std::vector<std::string> iso_wp_list = {
        ::cuts::electronID_Run2::mvaEleID_iso_Tight, ::cuts::electronID_Run2::mvaEleID_iso_Medium
    };
    for(const std::string& wp_name : iso_wp_list) {
        const bool l1_pass = l1->electronID(wp_name) > 0.5f;
        const bool l2_pass = l2->electronID(wp_name) > 0.5f;
        if(l1_pass && !l2_pass) return true;
        if(!l1_pass && l2_pass) return false;
        if(l1_pass && l2_pass) break;
    }
    return l1.IsMoreIsolated(l2);
}

}

class TupleStore {
public:
    using Mutex = std::recursive_mutex;

    static ntuple::EventTuple& GetTuple();
    static void ReleaseEventTuple();

private:
    static Mutex mutex;
    static std::atomic<int> tuple_counter;
    static std::shared_ptr<ntuple::EventTuple> eventTuple_ptr;
};

class BaseTupleProducer : public edm::EDAnalyzer {
public:
    using ElectronCandidate = analysis::SelectionResultsBase::ElectronCandidate;
    using MuonCandidate = analysis::SelectionResultsBase::MuonCandidate;
    using TauCandidate = analysis::SelectionResultsBase::TauCandidate;
    using JetCandidate = analysis::SelectionResultsBase::JetCandidate;
    using MET = analysis::MissingET<pat::MET>;
    using MetCovMatrix = MET::CovMatrix;
    using SelectionManager = analysis::SelectionManager;
    using Cutter = cuts::Cutter<SelectionManager>;
    using LorentzVector = analysis::LorentzVector;
    using LorentzVectorM = analysis::LorentzVectorM;
    using LorentzVectorE = analysis::LorentzVectorE;

    static constexpr double pt_shift = 5;

private:
    std::string treeName;
    TupleProducerData anaData;
    std::map<std::string, std::shared_ptr<SelectionData>> anaDataBeforeCut, anaDataAfterCut;
    edm::EDGetTokenT<std::vector<pat::Electron>> electronsMiniAOD_token;
    edm::EDGetTokenT<std::vector<pat::Tau>> tausMiniAOD_token;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsMiniAOD_token;
    edm::EDGetTokenT<edm::View<reco::Vertex>> vtxMiniAOD_token;
    edm::EDGetTokenT<edm::View<pat::MET>> pfMETAOD_token;
    edm::EDGetTokenT<edm::View<reco::GenMET>> genMETAOD_token;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsMiniAOD_token;
    edm::EDGetTokenT<std::vector<pat::Jet>> fatJetsMiniAOD_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUInfo_token;
    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_token;
    edm::EDGetTokenT<GenEventInfoProduct> genWeights_token;
    edm::EDGetTokenT<TtGenEvent> topGenEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<edm::View<reco::GenJet>> genJets_token;
    edm::EDGetTokenT<double> prefweight_token;
    edm::EDGetTokenT<double> prefweightup_token;
    edm::EDGetTokenT<double> prefweightdown_token;
    edm::EDGetTokenT<double> rho_token;
    std::map<std::string, edm::EDGetTokenT<bool>> customMetFilters_token;
    edm::EDGetTokenT<edm::ValueMap<float>> updatedPileupJetIdDiscr_token;
    edm::EDGetTokenT<edm::ValueMap<int>> updatedPileupJetId_token;

protected:
    const analysis::Period period;
    const bool isMC, applyTriggerMatch, applyTriggerMatchCut, applyTriggerCut, storeLHEinfo;
    const int nJetsRecoilCorr;
    const bool saveGenTopInfo, saveGenBosonInfo, saveGenJetInfo, saveGenParticleInfo, isEmbedded;
    ntuple::EventTuple& eventTuple;
    analysis::TriggerTools triggerTools;

private:
    const edm::Event *edmEvent;
    edm::Handle<std::vector<pat::Electron> > pat_electrons;
    edm::Handle<std::vector<pat::Tau> > pat_taus;
    edm::Handle<std::vector<pat::Muon> > pat_muons;
    edm::Handle<edm::View<reco::Vertex> > vertices;
    edm::Handle<edm::View<pat::MET> > pfMETs;
    edm::Handle<std::vector<pat::Jet> > pat_jets;
    edm::Handle<std::vector<pat::Jet> > pat_fatJets;
    edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
    edm::Handle<LHEEventProduct> lheEventProduct;
    edm::Handle<edm::View<reco::GenMET>> genMET;
    edm::Handle<GenEventInfoProduct> genEvt;
    edm::Handle<TtGenEvent> topGenEvent;
    edm::Handle<double> rho;

    JME::JetResolution resolution;

protected:
    edm::EventID eventId;
    std::vector<ElectronCandidate> electrons;
    std::vector<MuonCandidate> muons;
    std::vector<TauCandidate> taus;
    std::vector<JetCandidate> jets, fatJets;
    std::shared_ptr<MET> met;
    edm::Ptr<reco::Vertex> primaryVertex;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    edm::Handle<edm::View<reco::GenJet>> genJets;
    edm::Handle<edm::ValueMap<float>> updatedPileupJetIdDiscr;
    edm::Handle<edm::ValueMap<int>> updatedPileupJetId;

private:
    void InitializeAODCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void InitializeCandidateCollections();
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    virtual void endJob() override;
    virtual void ProcessEvent(Cutter& cut) = 0;
    static double Isolation(const pat::Electron& electron);
    static double Isolation(const pat::Muon& muon);

    template<typename Input>
    void ConsumeTag(edm::EDGetTokenT<Input>& token, const edm::ParameterSet& cfg, const std::string& name)
    {
        token = mayConsume<Input>(cfg.getParameter<edm::InputTag>(name));
    }

public:
    BaseTupleProducer(const edm::ParameterSet& iConfig, analysis::Channel _channel);

protected:
    TupleProducerData& GetAnaData() { return anaData; }

    static bool PassPFTightId(const pat::Jet& pat_jet, analysis::Period period);

    void ApplyBaseSelection(analysis::SelectionResultsBase& selection);
    void FillEventTuple(const analysis::SelectionResultsBase& selection);
    void FillElectron(const analysis::SelectionResultsBase& selection);
    void FillMuon(const analysis::SelectionResultsBase& selection);
    void FillTau(const analysis::SelectionResultsBase& selection);
    void FillHiggsDaughtersIndexes(const analysis::SelectionResultsBase& selection, size_t shift);
    void FillLheInfo();
    void FillGenParticleInfo();
    void FillGenJetInfo();
    void FillLegGenMatch(const analysis::LorentzVectorXYZ& p4);
    void FillMetFilters(analysis::Period period);
    void FillOtherLeptons(const std::vector<ElectronCandidate>& other_electrons, const std::vector<MuonCandidate>& other_muons);

    std::vector<ElectronCandidate> CollectVetoElectrons(bool isTightSelection,
            const std::vector<const ElectronCandidate*>& signalElectrons);
    std::vector<MuonCandidate> CollectVetoMuons(bool isTightSelection,
            const std::vector<const MuonCandidate*>& signalMuons);

    std::vector<JetCandidate> CollectJets();

    void SelectVetoElectron(const ElectronCandidate& electron, Cutter& cut,
                            const std::vector<const ElectronCandidate*>& signalElectrons,
                            bool isTightSelection) const;
    void SelectVetoMuon(const MuonCandidate& muon, Cutter& cut,
                        const std::vector<const MuonCandidate*>& signalMuons,
                        bool isTightSelection) const;

    void SelectJet(const JetCandidate& jet, Cutter& cut) const;
    bool PassMatchSelection(const TauCandidate& tau) const;
    bool PassIsoSelection(const TauCandidate& tau) const;

    template<typename Candidate1, typename Candidate2>
    std::vector<std::pair<size_t,size_t>> FindCompatibleObjects(const std::vector<Candidate1>& objects1,
                                                       const std::vector<Candidate2>& objects2,
                                                       double minDeltaR, const std::string& hist_name,
                                                       int expectedCharge = analysis::AnalysisObject::UnknownCharge)
    {
        const double minDeltaR2 = std::pow(minDeltaR, 2);
        std::vector<std::pair<size_t,size_t>> higgses_indexes;
        for(size_t n = 0; n < objects1.size(); ++n) {
            for(size_t m = 0; m < objects2.size(); ++m) {
                const auto& object1 = objects1.at(n);
                const auto& object2 = objects2.at(m);
                if(ROOT::Math::VectorUtil::DeltaR2(object1.GetMomentum(), object2.GetMomentum()) > minDeltaR2) {
                    const analysis::CompositeCandidate<Candidate1, Candidate2> candidate(object1, object2);
                    if (expectedCharge !=analysis::AnalysisObject::UnknownCharge
                            && candidate.GetCharge() != expectedCharge) continue;
                    higgses_indexes.emplace_back(n,m);
                    GetAnaData().Mass(hist_name).Fill(candidate.GetMomentum().M(), 1);
                }
            }
        }

        GetAnaData().N_objects(hist_name).Fill(higgses_indexes.size(), 1);
        return higgses_indexes;
    }

    template<typename Candidate, typename BaseSelectorType, typename Comparitor = std::less<Candidate>>
    std::vector<Candidate> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                              const std::vector<Candidate>& all_candidates,
                                              Comparitor comparitor = Comparitor())
    {
        static constexpr double weight = 1;
        const std::string suffix = selection_label;
        auto& objectSelector = GetAnaData().Selection(suffix);
        objectSelector.SetSave(true);

        if(!anaDataBeforeCut.count(suffix) )
              anaDataBeforeCut[suffix] = std::make_shared<SelectionData>(&edm::Service<TFileService>()->file(),
                                                                         treeName + "_before_cut/" + suffix);

        auto entry = &anaDataBeforeCut.at(suffix)->h;
        SelectionManager selectionManager(entry, weight, selection_label);

        const auto selector = [&](size_t id) -> Candidate {
            const Candidate& candidate = all_candidates.at(id);

            Cutter cut(&objectSelector, &selectionManager);
            base_selector(candidate, cut);
            return candidate;
        };

        const auto selected = objectSelector.template collect_objects<Candidate>(1, all_candidates.size(), selector,
                                                                                 comparitor);

        if(!anaDataAfterCut.count(suffix))
            anaDataAfterCut[suffix] = std::make_shared<SelectionData>(&edm::Service<TFileService>()->file(),
                                                                       treeName + "_after_cut/" + suffix);

        auto entry_afterCut = &anaDataAfterCut.at(suffix)->h;
        SelectionManager selectionManager_afterCut(entry_afterCut, weight, selection_label);

        for(const auto& candidate : selected) {
            Cutter cut(nullptr, &selectionManager_afterCut);
            base_selector(candidate, cut);
        }

        GetAnaData().N_objects(suffix).Fill(selected.size(), weight);

        return selected;
    }

    template<typename Candidate>
    static float GetUserFloat(const Candidate& obj, const std::string& name)
    {
        return obj->userFloat(name);
    }
};
