/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define SUMMARY_DATA() \
    VAR(ULong64_t, numberOfProcessedEvents) \
    VAR(Double_t, totalWeight) \
    VAR(std::vector<std::string>, tauId_names) \
    VAR(std::vector<uint32_t>, tauId_keys) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ProdSummary, SummaryTuple, SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, SummaryTuple, SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA
