/*! Definition of the h->tautau sync tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define TTBar_DATA() \
    VAR(Int_t, run) /* Run */ \
    VAR(Int_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    VAR(Int_t, channelID) /* Channel: MuTau, ETau, TauTau */ \
    VAR(Int_t, eventType) /* event type category */ \
    /*Top variables*/ \
    VAR(std::vector<Int_t>, pdgid_tops) /* pdgid  top*/ \
    VAR(std::vector<Float_t>, pt_tops) /* pt  top*/ \
    VAR(std::vector<Float_t>, eta_tops) /* eta  top*/ \
    VAR(std::vector<Float_t>, phi_tops) /* phi  top*/ \
    VAR(std::vector<Float_t>, mass_tops) /* mass  top*/ \
    VAR(std::vector<Float_t>, energy_tops) /* energy  top*/ \
    /* W variables */  \
    VAR(std::vector<Int_t>, DecayMode_Ws) /* pdgid  W*/ \
    VAR(std::vector<Int_t>, pdgid_Ws) /* pdgid  W*/ \
    VAR(std::vector<Float_t>, pt_Ws) /* pt  W*/ \
    VAR(std::vector<Float_t>, eta_Ws) /* eta  W*/ \
    VAR(std::vector<Float_t>, phi_Ws) /* phi  W*/ \
    VAR(std::vector<Float_t>, mass_Ws) /* mass  W*/ \
    VAR(std::vector<Float_t>, energy_Ws) /* energy  W*/ \
    /* First Ws daugther*/ \
    VAR(std::vector<Int_t>, pdgid_dau1) /* pdgid  */ \
    VAR(std::vector<Float_t>, pt_dau1) /* pt  */ \
    VAR(std::vector<Float_t>, eta_dau1) /* eta  */ \
    VAR(std::vector<Float_t>, phi_dau1) /* phi  */ \
    VAR(std::vector<Float_t>, mass_dau1) /* mass  */ \
    VAR(std::vector<Float_t>, energy_dau1) /* energy  */ \
    /* Second Ws daugther*/ \
    VAR(std::vector<Int_t>, pdgid_dau2) /* pdgid  */ \
    VAR(std::vector<Float_t>, pt_dau2) /* pt  */ \
    VAR(std::vector<Float_t>, eta_dau2) /* eta  */ \
    VAR(std::vector<Float_t>, phi_dau2) /* phi  */ \
    VAR(std::vector<Float_t>, mass_dau2) /* mass  */ \
    VAR(std::vector<Float_t>, energy_dau2) /* energy  */ \
    /* BJets */ \
    VAR(std::vector<Int_t>, pdgid_bjet) /* pdgid  */ \
    VAR(std::vector<Float_t>, pt_bjet) /* pt  */ \
    VAR(std::vector<Float_t>, eta_bjet) /* eta  */ \
    VAR(std::vector<Float_t>, phi_bjet) /* phi  */ \
    VAR(std::vector<Float_t>, mass_bjet) /* mass  */ \
    VAR(std::vector<Float_t>, energy_bjet) /* energy  */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, TTBar, TTBarTree, TTBar_DATA, "ttbarStudy")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, TTBarTree, TTBar_DATA)
#undef VAR
#undef SYNC_DATA

namespace Run2 {
inline double DefaultFillValueForSyncTree() { return -10000; }
inline float DefaultFloatFillValueForSyncTree() { return std::numeric_limits<float>::lowest(); }

enum class WDecayMode { Quarks = 0, Mu = 1, Ele = 2, Tau = 3};

}
