/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define SUMMARY_DATA() \
    /* Run statistics */ \
    VAR(UInt_t, exeTime) \
    VAR(ULong64_t, numberOfProcessedEvents) \
    VAR(Double_t, totalWeight) \
    /* Tau ID information */ \
    VAR(std::vector<std::string>, tauId_names) \
    VAR(std::vector<uint32_t>, tauId_keys) \
    /* MC truth event splitting */ \
    VAR(std::vector<UInt_t>, lhe_n_partons) \
    VAR(std::vector<UInt_t>, lhe_n_b_partons) \
    VAR(std::vector<UInt_t>, lhe_ht10_bin) \
    VAR(std::vector<ULong64_t>, lhe_n_events) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ProdSummary, SummaryTuple, SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, SummaryTuple, SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA

namespace ntuple {
struct GenId {
    size_t n_partons;
    size_t n_b_partons;
    size_t ht10_bin;

    GenId() : n_partons(0), n_b_partons(0), ht10_bin(0) {}
    GenId(size_t _n_partons, size_t _n_b_partons, size_t _ht10_bin) :
        n_partons(_n_partons), n_b_partons(_n_b_partons), ht10_bin(_ht10_bin) {}

    bool operator<(const GenId& other) const
    {
        if(n_partons != other.n_partons) return n_partons < other.n_partons;
        if(n_b_partons != other.n_b_partons) return n_b_partons < other.n_b_partons;
        return ht10_bin < other.ht10_bin;
    }
};

using GenEventCountMap = std::map<GenId, size_t>;

} // namespace ntuple
