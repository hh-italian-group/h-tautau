/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "EventTuple.h"

namespace cache_tuple {
using LorentzVectorM = ntuple::LorentzVectorM;
}

#define CACHE_DATA() \
    VAR(UInt_t, run) /* run */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Long64_t, entry_index) /* entry index in the original tuple */ \
    /* SVfit variables */ \
    VAR(std::vector<UInt_t>, SVfit_htt_index) /* SVfit: Higgs index */ \
    VAR(std::vector<Bool_t>, SVfit_is_valid) /* SVfit: has a valid result */ \
    VAR(std::vector<LorentzVectorM>, SVfit_p4) /* SVfit: 4-momentum */ \
    VAR(std::vector<LorentzVectorM>, SVfit_p4_error) /* SVfit: error on 4-momentum */ \
    VAR(std::vector<Float_t>, SVfit_mt) /* SVfit: transverse mass */ \
    VAR(std::vector<Float_t>, SVfit_mt_error) /* SVfit: error on transverse mass */ \
    VAR(std::vector<Int_t>, SVfit_unc_source) /* SVfit: uncertainty source */ \
    VAR(std::vector<Int_t>, SVfit_unc_scale) /* SVfit: uncertainty scale */ \
    /* HHKinFit variables */ \
    VAR(std::vector<UInt_t>, kinFit_htt_index) /* HHKinFit: H->tautau index */ \
    VAR(std::vector<UInt_t>, kinFit_hbb_index) /* HHKinFit: H->bb index */\
    VAR(std::vector<Int_t>, kinFit_unc_source) /* HHKinFit: uncertianty source */ \
    VAR(std::vector<Int_t>, kinFit_unc_scale) /* HHKinFit: uncertainty scale */ \
    VAR(std::vector<Float_t>, kinFit_m) /* HHKinFit: m_bbtt mass */\
    VAR(std::vector<Float_t>, kinFit_chi2) /*  HHKinFit: chi2 value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* HHKinFit: convergence code */\
    /* Jet HH-btag score variables */ \
    VAR(std::vector<UInt_t>, jet_HHbtag_htt_index) /* HH-btag: H->tautau index */ \
    VAR(std::vector<UInt_t>, jet_HHbtag_jet_index) /* HH-btag: jet index */ \
    VAR(std::vector<Int_t>, jet_HHbtag_unc_source) /* HH-btag: uncertainty source */ \
    VAR(std::vector<Int_t>, jet_HHbtag_unc_scale) /* HH-btag: uncertainty scale */ \
    VAR(std::vector<Float_t>, jet_HHbtag_value) /* HH-btag: tagging score */ \
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
DECLARE_TREE(cache_tuple, CacheProdSummary, CacheSummaryTuple, CACHE_SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(cache_tuple, CacheSummaryTuple, CACHE_SUMMARY_DATA)
#undef VAR
#undef CACHE_SUMMARY_DATA
