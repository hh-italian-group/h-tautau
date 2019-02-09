/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "EventTuple.h"
#include "AnalysisTypes.h"

#define SUMMARY_DATA() \
    /* Run statistics */ \
    VAR(UInt_t, exeTime) \
    VAR(ULong64_t, numberOfProcessedEvents) \
    VAR(Double_t, totalShapeWeight) \
    VAR(Double_t, totalShapeWeight_withTopPt) \
    /* Tau ID information */ \
    VAR(std::vector<std::string>, tauId_names) \
    /* Trigger information */ \
    VAR(std::vector<Int_t>, triggers_channel) \
    VAR(std::vector<std::string>, triggers_pattern) \
    /* Skimmer Variables */\
    VAR(std::vector<UInt_t>, file_desc_id) /* vector of File id in TupleSkimmer. */ \
    VAR(std::vector<std::string>, file_desc_name) /* vector of File name in TupleSkimmer. */ \
    VAR(UInt_t, n_splits) /* Number of splits for a file in TupleSkimmer. */ \
    VAR(UInt_t, split_seed) /* Seed for splitting in TupleSkimmer. */ \
    /* MC truth event splitting */ \
    VAR(std::vector<UInt_t>, lhe_n_partons) \
    VAR(std::vector<UInt_t>, lhe_n_b_partons) \
    VAR(std::vector<UInt_t>, lhe_ht10_bin) \
    VAR(std::vector<ULong64_t>, lhe_n_events) \
    /* Top reweighting */ \
    VAR(std::vector<Int_t>, genEventType) /* top gen event type */ \
    VAR(std::vector<ULong64_t>, genEventType_n_events) /* n events for top gen event type */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ProdSummary, SummaryTuple, SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, SummaryTuple, SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA


#define EVENT_EXPRESS_DATA() \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(UInt_t, run) /* run */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, gen_top_pt) /* pt of gen ME top */ \
    VAR(Float_t, gen_topBar_pt) /* pt of gen ME anti-top */ \
    VAR(Int_t,   genEventType) /* top gen event type */ \
    VAR(Float_t, lhe_H_m) /* mass of lhe H */ \
    VAR(Float_t, lhe_hh_m) /* mass of lhe hh pair */ \
    VAR(Float_t, lhe_hh_cosTheta) /* cos(theta) between h and z-axis in the hh reference frame */ \
    /* MC truth event splitting */ \
    VAR(UInt_t, file_desc_id)  /*File id in ExpressTuple.*/  \
    VAR(UInt_t, split_id) /* Split id in TupleSkimmer. */ \
    VAR(UInt_t, lhe_n_partons) \
    VAR(UInt_t, lhe_n_b_partons) \
    VAR(UInt_t, lhe_ht10_bin) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ExpressEvent, ExpressTuple, EVENT_EXPRESS_DATA, "all_events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, ExpressTuple, EVENT_EXPRESS_DATA)
#undef VAR
#undef SUMMARY_DATA

namespace ntuple {
struct GenId {
    size_t n_partons;
    size_t n_b_partons;
    size_t ht10_bin;

    GenId();
    GenId(size_t _n_partons, size_t _n_b_partons, size_t _ht10_bin);
    bool operator<(const GenId& other) const;
};

using GenEventCountMap = std::map<GenId, size_t>;
using GenEventTypeCountMap = std::map<analysis::GenEventType, size_t>;

GenEventCountMap ExtractGenEventCountMap(const ProdSummary& s);
void ConvertGenEventCountMap(ProdSummary& s, const GenEventCountMap& genCountMap);
GenEventTypeCountMap ExtractGenEventTypeCountMap(const ProdSummary& s);
void ConvertGenEventTypeCountMap(ProdSummary& s, const GenEventTypeCountMap& genCountMap);
std::shared_ptr<SummaryTuple> CreateSummaryTuple(const std::string& name, TDirectory* directory,
                                                 bool readMode, TreeState treeState);
void CheckProdSummaryConsistency(const ProdSummary& s);
bool CheckProdSummaryCompatibility(const ProdSummary& s1, const ProdSummary& s2, std::ostream* os = nullptr);
void MergeProdSummaries(ProdSummary& summary, const ProdSummary& otherSummary);
ProdSummary MergeSummaryTuple(SummaryTuple& tuple);
std::shared_ptr<ExpressTuple> CreateExpressTuple(const std::string& name, TDirectory* directory,
                                                 bool readMode, TreeState treeState);
} // namespace ntuple
