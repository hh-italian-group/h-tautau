/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define CACHE_SUMMARY_DATA() \
    /* Run statistics */ \
    VAR(UInt_t, exeTime) \
    VAR(Int_t, numberOfOriginalEvents) \
    VAR(Int_t, numberOfTimesSVFit) \
    VAR(Int_t, numberOfTimesKinFit) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(cache_ntuple, CacheProdSummary, CacheSummaryTuple, CACHE_SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(cache_ntuple, CacheSummaryTuple, CACHE_SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA
