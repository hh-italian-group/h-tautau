/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

#define LVAR(type, name, n) VAR(type, name##_##n)

#define LEG_DATA(n) \
    LVAR(analysis::LorentzVectorM, p4, n) /* 4-momentum */ \
    LVAR(Int_t, q, n) /* Charge */ \
    LVAR(Float_t, d0, n) /* d0 with respect to primary vertex */ \
    LVAR(Float_t, dZ, n) /* dZ with respect to primary vertex */ \
    LVAR(Float_t, iso, n) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    LVAR(Float_t, id_e_mva_nt_loose, n) /* Non-triggering electron ID MVA score id (when using electron) 0 otherwise */ \
    LVAR(Int_t, gen_match, n) /*Generator matching, see Htautau Twiki*/\
    LVAR(root_ext::strmap<float>, tauIDs, n) /* tau ID variables */ \
    /**/

#define JVAR(type, name, col) VAR(std::vector<type>, col##_##name)

#define JET_COMMON_DATA(col) \
    JVAR(analysis::LorentzVectorE, p4, col) /* Jet 4-momentum */ \
    JVAR(Float_t, csv, col) /* Jet CSV value */ \
    /**/

#define JET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(Float_t, rawf, col) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Float_t, mva, col) /* Jet MVA id value */ \
    JVAR(Int_t, hadronFlavour, col) \
    /**/

#define FATJET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(Float_t, m_pruned, col) \
    JVAR(Float_t, m_filtered, col) \
    JVAR(Float_t, m_trimmed, col) \
    JVAR(Float_t, m_softDrop, col) \
    JVAR(Float_t, n_subjettiness_tau1, col) \
    JVAR(Float_t, n_subjettiness_tau2, col) \
    JVAR(Float_t, n_subjettiness_tau3, col) \
    /**/

#define SUBJET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(size_t, parentIndex, col) \
    /**/

#define MVAR(type, name, col) VAR(type, col##_##name)

#define MET_DATA(col) \
    MVAR(analysis::LorentzVectorM, p4, col) /* MET 4-momentum */ \
    MVAR(analysis::SquareMatrix<2>, cov, col) /* pf met covariance matrix */ \
    /**/

#define EVENT_DATA() \
    VAR(UInt_t, run) /* Run */ \
    VAR(UInt_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    VAR(Int_t, channelID) /* Channel: MuTau, ETau, TauTau */ \
    VAR(Int_t, eventEnergyScale) /* event type category */ \
    VAR(Int_t, eventType) /* event type category */ \
    VAR(Double_t, weightevt) /*Gen Event Weight*/ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* NPV */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    /* SV Fit variables */ \
    VAR(analysis::LorentzVectorM, SVfit_p4) /* SV Fit using integration method */ \
    /* Signal leptons */ \
    LEG_DATA(1) /* muon for muTau, electron for eTau, electron for eMu, Leading (in pT) Tau for tauTau */ \
    LEG_DATA(2) /* hadronic Tau for muTau and eTau, Muon for eMu, Trailing (in pT) Tau for tauTau */ \
    /* Met related variables */ \
    MET_DATA(pfMET) \
    MET_DATA(mvaMET) \
    MET_DATA(puppiMET) \
    /* Candidate Jets: jets after applying Jet energy corrections (excluding hadronic Tau) */ \
    JET_DATA(jets) \
    FATJET_DATA(fatJets) \
    SUBJET_DATA(subJets) \
    /* KinFit Variables */ \
    VAR(std::vector<Double_t>, kinFit_m) /* KinFit m_bbtt mass compute the first 2 jets, ordered by CSV*/\
    VAR(std::vector<Double_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Double_t>, kinFit_probability) /*  KinFit chi2 probability value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    /* LHE info */\
    VAR(std::vector<Int_t>, lhe_particle_pdg) \
    VAR(std::vector<analysis::LorentzVectorM>, lhe_particle_p4) \
    VAR(UInt_t, lhe_n_partons) \
    VAR(UInt_t, lhe_n_b_partons) \
    VAR(Float_t, lhe_HT) \
    /* Vetos */\
    VAR(Bool_t, dilepton_veto) /* Event is vetoed by the dilepton veto if true */ \
    VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, Event, EventTuple, EVENT_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, EventTuple, EVENT_DATA)
#undef VAR
#undef EVENT_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef FATJET_DATA
#undef SUBJET_DATA
#undef JET_COMMON_DATA
#undef JVAR
#undef MET_DATA
#undef MVAR

namespace ntuple {
template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }
}
