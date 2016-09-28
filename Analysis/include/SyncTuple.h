/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

#define LVAR(type, name, pref) VAR(type, pref##_##name)

#define LEG_DATA(pref) \
        LVAR(Float_t, iso, pref) /* iso */ \
        LVAR(Int_t, MVAiso, pref) /* mva Iso working points */ \
        LVAR(Float_t, pt, pref) \
        LVAR(Float_t, eta, pref) \
        LVAR(Float_t, phi, pref) \
        LVAR(Float_t, e, pref) \
        LVAR(Float_t, mT, pref) \

#define JET_DATA(pref) \
        LVAR(Float_t, pt, pref) \
        LVAR(Float_t, eta, pref) \
        LVAR(Float_t, phi, pref) \
        LVAR(Float_t, e, pref) \
        LVAR(Int_t, bID, pref) \
        LVAR(Int_t, flav, pref) \
        LVAR(Float_t, pt_raw, pref) \


#define SYNC_DATA() \
        VAR(Int_t, RunNumber) /* Run */ \
        VAR(Int_t, lumi) /* Lumi */ \
        VAR(ULong64_t, EventNumber) /* Evt */ \
        VAR(Int_t, pairType) /* Channel: MuTau, ETau, TauTau */ \
        VAR(Int_t, eventEnergyScale) /* event type category */ \
        /* Event Variables */ \
        VAR(Int_t, npv) /* NPV */ \
        VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
        VAR(Float_t, met_phi)  \
        VAR(Float_t, met_et)  \
        VAR(Float_t, met_et_corr)  \
        LEG_DATA(dau1) \
        LEG_DATA(dau2) \
        JET_DATA(bjet1) \
        JET_DATA(bjet2) \

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, SyncEvent, SyncTuple, SYNC_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, SyncTuple, SYNC_DATA)
#undef VAR
#undef SYNC_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef JVAR
