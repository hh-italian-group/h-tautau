/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace cache_tuple {
using LorentzVectorM = analysis::LorentzVectorM;
}

#define CACHE_DATA() \
    VAR(UInt_t, run) /* run */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    /* SV Fit variables */ \
    VAR(std::vector<size_t>, SVfit_Higgs_index) /* SVfit using integration method */ \
    VAR(std::vector<Bool_t>, SVfit_is_valid) /* SVfit using integration method */ \
    VAR(std::vector<LorentzVectorM>, SVfit_p4) /* SVfit using integration method */ \
    VAR(std::vector<LorentzVectorM>, SVfit_p4_error) /* SVfit using integration method */ \
    VAR(std::vector<Float_t>, SVfit_mt) /* SVfit using integration method */ \
    VAR(std::vector<Float_t>, SVfit_mt_error) /* SVfit using integration method */ \
    VAR(std::vector<Int_t>, SVfit_unc_source) /* SVfit using integration method */ \
    VAR(std::vector<Int_t>, SVfit_unc_scale) /* SVfit using integration method */ \
    /* KinFit Variables */ \
    VAR(std::vector<size_t>, kinFit_Higgs_index) /* kinFit Higgs indexes */ \
    VAR(std::vector<size_t>, kinFit_jetPairId) /* indices of jet pairs for which KinFit is calculated */\
    VAR(std::vector<Float_t>, kinFit_m) /* KinFit m_bbtt mass */\
    VAR(std::vector<Float_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    VAR(std::vector<Int_t>, kinFit_unc_source) /* kinFit */ \
    VAR(std::vector<Int_t>, kinFit_unc_scale) /* kinFit */ \
    /* Jet score Variables */ \
    VAR(std::vector<size_t>, jet_hh_score_index) /* jet score indexes */ \
    VAR(std::vector<Int_t>, jet_hh_score_unc_scale) /* jet score scale */ \
    VAR(std::vector<Int_t>, jet_hh_score_unc_source) /* jet score source */ \
    VAR(std::vector<Float_t>, jet_hh_score_value) /* jet score value */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(cache_tuple, CacheEvent, CacheTuple, CACHE_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(cache_tuple, CacheTuple, CACHE_DATA)
#undef VAR
#undef CACHE_DATA


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
#undef CACHE_SUMMARY_DATA
