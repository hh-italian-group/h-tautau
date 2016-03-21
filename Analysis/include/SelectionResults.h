/*! Definition of selection results container used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "Candidate.h"

#include "SVfit.h"
//#include "KinFit.h"

#define SELECTION_ENTRY(name) \
    ANA_DATA_ENTRY(cuts::ObjectSelector, name) \
    /**/

#define X(name) \
    selectionManager.FillHistogram(object.name, #name)

#define Y(name) \
    selectionManager.FillHistogram(name, #name)


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

//struct SelectionResults {
//    static constexpr size_t NumberOfLegs = 2;

//    virtual ~SelectionResults() {}
//    CandidatePtr higgs;
//    sv_fit::CombinedFitResults svfitResults;
//    kinematic_fit::four_body::FitResults kinfitResults;
//    CandidatePtrVector jets;
//    CandidatePtrVector jetsPt20;
//    CandidatePtrVector bjets_all;
//    CandidatePtrVector retagged_bjets;
//    VertexPtrVector vertices;
//    ntuple::MET pfMET;
//    ntuple::MET MET_with_recoil_corrections;
//    ntuple::EventType eventType;

//    VertexPtr GetPrimaryVertex() const { return vertices.front(); }
//    virtual CandidatePtr GetLeg(size_t leg_id) const = 0;
//    virtual const finalState::bbTauTau& GetFinalStateMC() const = 0;
//};

struct SelectionResultsV2 {
    static constexpr size_t NumberOfLegs = 2;

    virtual ~SelectionResultsV2() {}
    CandidateV2Ptr higgs;
    Float_t numtruepileupinteractions =-1;
    Float_t HT;
    Double_t weightevt;
    bool Zveto;
    sv_fit::FitResults svfitResult;
    //kinematic_fit::four_body::FitResults kinfitResults;
    CandidateV2PtrVector jets;
    int numJet;
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
} // namespace analysis
