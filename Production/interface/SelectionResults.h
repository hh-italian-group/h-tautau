/*! Definition of selection results container used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/CutTools.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/Candidate.h"
#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "SVfitInterface.h"

#define SELECTION_ENTRY(name) \
    ANA_DATA_ENTRY(cuts::ObjectSelector, name) \
    /**/

namespace analysis {
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

struct SelectionResultsBase {
    static constexpr size_t NumberOfLegs = 2;
    using TauCandidate = LeptonCandidate<pat::Tau>;
    using JetCandidate = Candidate<pat::Jet>;
    using JetCandidateVector = std::vector<JetCandidate>;
    using Vertex = reco::Vertex;

    bool Zveto, electronVeto, muonVeto;
    sv_fit::FitResults svfitResult;
    std::vector<kin_fit::FitResults> kinfitResults;
    JetCandidateVector jets;
    JetCandidateVector bjets;
    const Vertex* primaryVertex;

    virtual const TauCandidate& GetSecondLeg() const = 0;
    virtual const LorentzVector& GetHiggsMomentum() const = 0;
    virtual ~SelectionResultsBase() {}
};

template<typename _FirstLeg>
struct SelectionResults : SelectionResultsBase {
    using FirstLeg = _FirstLeg;
    using SecondLeg = TauCandidate;
    using HiggsCandidate = CompositCandidate<FirstLeg, SecondLeg>;
    using HiggsCandidatePtr = std::shared_ptr<HiggsCandidate>;

    HiggsCandidatePtr higgs;

    void SetHiggsCandidate(const HiggsCandidate& h) { higgs = HiggsCandidatePtr(new HiggsCandidate(h)); }
    virtual const TauCandidate& GetSecondLeg() const override { return higgs->GetSecondDaughter(); }
    virtual const LorentzVector& GetHiggsMomentum() const override { return higgs->GetMomentum(); }
};

} // namespace analysis
