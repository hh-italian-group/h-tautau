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
    LVAR(Float_t, mt, n) /* mT of the lepton wrt to MVA met */ \
    LVAR(Float_t, pfmt, n) /* mT of  first lepton wrt to PF met */ \
    LVAR(Float_t, puppimt, n) /* mT of  first lepton wrt to PUPPI met */ \
    LVAR(Float_t, iso, n) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    LVAR(Float_t, id_e_mva_nt_loose, n) /* Non-triggering electron ID MVA score id (when using electron) 0 otherwise */ \
    LVAR(Int_t, gen_match, n) /*Generator matching, see Htautau Twiki*/\
    LVAR(Float_t, againstElectronLooseMVA6, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstElectronMediumMVA6, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstElectronTightMVA6, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstElectronVLooseMVA6, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstElectronVTightMVA6, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, byIsolationMVA2raw, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstMuonLoose3, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, againstMuonTight3, n) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    LVAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits, n) /*  */ \
    LVAR(Float_t, byIsolationMVA3newDMwoLTraw, n) /* MVA iso for Tau w/o Lifetime information New Decay Mode */ \
    LVAR(Float_t, byIsolationMVA3oldDMwoLTraw, n) /* MVA iso for Tau w/o Lifetime information Old Decay Mode */ \
    LVAR(Float_t, byIsolationMVA3newDMwLTraw, n) /* MVA iso for Tau w/ Lifetime information New Decay Mode */ \
    LVAR(Float_t, byIsolationMVA3oldDMwLTraw, n) /* MVA iso for Tau w/ Lifetime information Old Decay Mode */ \
    LVAR(Int_t, byVLooseIsolationMVArun2v1DBoldDMwLT, n) /* MVA iso for Tau w/ Lifetime VLoose WP */ \
    LVAR(Int_t, byLooseIsolationMVArun2v1DBoldDMwLT, n) /* MVA iso for Tau w/ Lifetime Loose WP */ \
    LVAR(Int_t, byMediumIsolationMVArun2v1DBoldDMwLT, n) /* MVA iso for Tau w/ Lifetime Medium WP */ \
    LVAR(Int_t, byTightIsolationMVArun2v1DBoldDMwLT, n) /* MVA iso for Tau w/ Lifetime Tight WP */ \
    LVAR(Int_t, byVTightIsolationMVArun2v1DBoldDMwLT, n) /* MVA iso for Tau w/ Lifetime VTight WP */ \
    LVAR(Int_t, decayModeFindingOldDMs, n) /* Old Decay Mode finding */\
    /**/

#define JVAR(type, name, col) VAR(std::vector<type>, col##_##name)

#define JET_DATA(col) \
    JVAR(analysis::LorentzVectorE, p4, col) /* Jet 4-momentum */ \
    JVAR(Float_t, rawf, col) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Float_t, mva, col) /* Jet MVA id value */ \
    JVAR(Float_t, csv, col) /* Jet CSV value */ \
    JVAR(Int_t, partonFlavour, col) \
    /**/

#define MVAR(type, name, col) VAR(type, col##_##name)

#define MET_DATA(col) \
    MVAR(analysis::LorentzVectorM, p4, col) /* MET 4-momentum */ \
    MVAR(analysis::SquareMatrix<2>, cov, col) /* pf met covariance matrix */ \
    /**/

#define EVENT_DATA() \
    VAR(Int_t, run) /* Run */ \
    VAR(Int_t, lumi) /* Lumi */ \
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
    JET_DATA(fatjets) \
    /* KinFit Variables */ \
    VAR(std::vector<Double_t>, kinFit_m) /* KinFit m_bbtt mass compute the first 2 jets, ordered by CSV*/\
    VAR(std::vector<Double_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Double_t>, kinFit_probability) /*  KinFit chi2 probability value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    /* LHE info */\
    VAR(std::vector<Int_t>, lhe_particle_pdg) \
    VAR(std::vector<analysis::LorentzVectorM>, lhe_particle_p4) \
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
#undef JVAR
#undef MET_DATA
#undef MVAR


namespace ntuple {
template<typename T>
T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }

enum class HTbinning { lt0 = -1, lt100 = 0, f100to200 = 1, f200to400 = 2, f400to600 = 3, gt600 = 4 };
inline HTbinning GetHTbin(double HT)
{
    if(HT < 0) return HTbinning::lt0;
    if(HT < 100) return HTbinning::lt100;
    if(HT < 200) return HTbinning::f100to200;
    if(HT < 400) return HTbinning::f200to400;
    if(HT < 600) return HTbinning::f400to600;
    return HTbinning::gt600;
}

}
