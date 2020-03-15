/*! Definition of selection results container used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/CutTools.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/TriggerResults.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"

#define SELECTION_ENTRY(name) \
    ANA_DATA_ENTRY(cuts::ObjectSelector, name) \
    /**/

namespace analysis {
class SelectionManager {
public:
    using RootHist = TH1D;
    using Hist = root_ext::SmartHistogram<RootHist>;
    using Entry = root_ext::AnalyzerDataEntry<RootHist>;
    using AnaData = root_ext::AnalyzerData;
    using Factory = root_ext::HistogramFactory<RootHist>;
    SelectionManager(Entry* _entry, double _weight, const std::string& _selection_label)
        : entry(_entry), weight(_weight), selection_label(_selection_label) {}

    template<typename ValueType>
    ValueType FillHistogram(ValueType value, const std::string& hist_name)
    {
        if(entry != nullptr) {
            if(!entry->GetHistograms().count(hist_name)) {
                std::shared_ptr<Hist> hist(Factory::Make(hist_name, selection_label));
                entry->Set(hist_name, hist);
            }
            (*entry)(hist_name).Fill(value, weight);
        }
        return value;
    }

private:
    Entry* entry;
    double weight;
    std::string selection_label;
};

struct SelectionResultsBase {
    static constexpr size_t NumberOfLegs = 2;
    using TauCandidate = LeptonCandidate<pat::Tau>;
    using TauCandidateVector = std::vector<TauCandidate>;
    using JetCandidate = Candidate<pat::Jet, edm::Ptr<pat::Jet>>;
    using JetCandidateVector = std::vector<JetCandidate>;
    using ElectronCandidate = analysis::LeptonCandidate<pat::Electron, edm::Ptr<pat::Electron>>;
    using MuonCandidate = analysis::LeptonCandidate<pat::Muon>;
    using ElectronCandidateVector = std::vector<ElectronCandidate>;
    using MuonCandidateVector = std::vector<MuonCandidate>;
    using Vertex = reco::Vertex;

    edm::EventID eventId;

    bool Zveto, electronVeto, muonVeto;
    JetCandidateVector jets;
    TauCandidateVector taus;
    ElectronCandidateVector electrons;
    MuonCandidateVector muons;
    ElectronCandidateVector other_electrons;
    MuonCandidateVector other_muons;
    const Vertex* primaryVertex;
    std::vector<TriggerResults> triggerResults;
    std::vector<std::pair<size_t,size_t>> higgses_pair_indexes;

    SelectionResultsBase(const edm::EventID& _eventId) : eventId(_eventId) {}

    virtual ~SelectionResultsBase() {}
    // virtual const LorentzVector& GetHiggsMomentum() const = 0;

    bool HaveSameJets(const SelectionResultsBase& other) const
    {

        if(eventId != other.eventId) return false;
        if(jets.size() != other.jets.size()) return false;

        for(size_t n = 0; n < jets.size(); ++n) {
            if(&(*jets.at(n)) != &(*other.jets.at(n)))
                return false;
        }
        return true;
    }

    bool HaveSameOtherLeptons(const SelectionResultsBase& other) const
    {
        if(eventId != other.eventId) return false;
        if(other_electrons.size() != other.other_electrons.size()
                || other_muons.size() != other.other_muons.size()) return false;

        for(size_t n = 0; n < other_electrons.size(); ++n) {
            if(&(*other_electrons.at(n)) != &(*other.other_electrons.at(n)))
                return false;
        }

        for(size_t n = 0; n < other_muons.size(); ++n) {
            if(&(*other_muons.at(n)) != &(*other.other_muons.at(n)))
                return false;
        }
        return true;
    }
};

} // namespace analysis
