/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

using LorentzVectorM = analysis::LorentzVectorM_Float;

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
    VAR(std::vector<UInt_t>, kinFit_jetPairId) /* indices of jet pairs for which KinFit is calculated */\
    VAR(std::vector<Float_t>, kinFit_m) /* KinFit m_bbtt mass */\
    VAR(std::vector<Float_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    VAR(std::vector<Int_t>, kinFit_unc_source) /* kinFit */ \
    VAR(std::vector<Int_t>, kinFit_unc_scale) /* kinFit */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(cache_tuple, CacheEvent, CacheTuple, CACHE_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(cache_tuple, CacheTuple, CACHE_DATA)
#undef VAR
#undef CACHE_DATA
