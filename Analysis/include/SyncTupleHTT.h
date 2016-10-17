/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

#define LVAR(type, name, pref) VAR(type, name##_##pref)
#define JVAR(type, name, suff, pref) VAR(type, suff##name##_##pref)

#define LEG_DATA(pref) \
    LVAR(Float_t, pt, pref) \
    LVAR(Float_t, phi, pref) \
    LVAR(Float_t, eta, pref) \
    LVAR(Float_t, m, pref) \
    LVAR(Float_t, q, pref) \
    LVAR(Float_t, d0, pref) /* is dxy between leading track and first PV */ \
    LVAR(Float_t, dZ, pref) /* dZ between leading track and first PV */ \
    LVAR(Float_t, mt, pref) /* Use MVAMet sqrt(2 * l_pt * met_pt * (1 - cos( d_phi(l, met) )) */ \
    LVAR(Float_t, pfmt, pref) /* As above but using PF Met */ \
    LVAR(Float_t, puppimt, pref) /* As above but using Puppi Met */ \
    LVAR(Float_t, iso, pref) /* iso */ \
    LVAR(Float_t, id_e_mva_nt_loose, pref) /* Non-triggering electron ID MVA score */ \
    LVAR(Float_t, gen_match, pref) /* Type of gen particle matched to object (see HiggsToTauTauWorking2016#MC_Matching) */ \
    LVAR(Float_t, againstElectronLooseMVA6, pref) \
    LVAR(Float_t, againstElectronMediumMVA6, pref) \
    LVAR(Float_t, againstElectronTightMVA6, pref) \
    LVAR(Float_t, againstElectronVLooseMVA6, pref) \
    LVAR(Float_t, againstElectronVTightMVA6, pref) \
    LVAR(Float_t, againstMuonLoose3, pref) \
    LVAR(Float_t, againstMuonTight3, pref) \
    LVAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits, pref) \
    LVAR(Float_t, byIsolationMVA3newDMwoLTraw, pref) \
    LVAR(Float_t, byIsolationMVA3oldDMwoLTraw, pref) \
    LVAR(Float_t, byIsolationMVA3newDMwLTraw, pref) \
    LVAR(Float_t, byIsolationMVA3oldDMwLTraw, pref) \
    LVAR(Float_t, chargedIsoPtSum, pref) \
    LVAR(Float_t, decayModeFindingOldDMs, pref) \
    LVAR(Float_t, neutralIsoPtSum, pref) \
    LVAR(Float_t, puCorrPtSum, pref) \
    LVAR(Float_t, trigweight, pref) \
    LVAR(Float_t, idisoweight, pref) \
    /**/

#define JET_DATA(suff, pref) \
    JVAR(Float_t, pt, suff, pref) \
    JVAR(Float_t, eta, suff, pref) \
    JVAR(Float_t, phi, suff, pref) \
    JVAR(Float_t, rawf, suff, pref) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Int_t, mva, suff, pref) /* pu Jet id score */ \
    JVAR(Float_t, csv, suff, pref) \
    /**/

#define SYNC_DATA() \
    VAR(UInt_t, run) /* Run */ \
    VAR(UInt_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* Number of offline primary vertices */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(Float_t, rho) /* Use fixedGridRhoFastjetAll */ \
    /* Signal leptons */ \
    LEG_DATA(1) /* Leg 1 (leading tau for tt, electon for et,em muon for mt) */ \
    LEG_DATA(2) /* Leg 2 (trailing tau for tt, tau for et, mt, muon for em) */ \
    /* di-tau system */ \
    VAR(Float_t, pt_tt) /* like HIG-13-004 (p4_l1+p4_l2+p4_MET)_T, use PF met */  \
    VAR(Float_t, mt_tot) /* Use MVA MET (see HiggsToTauTauWorking2016#Synchronisation_Ntuple) */  \
    VAR(Float_t, m_vis) \
    VAR(Float_t, m_sv) /* using MarkovChain MC integration svfit.mass() */ \
    VAR(Float_t, mt_sv) /* using MarkovChain MC integration svfit.transverseMass() */ \
    /* MET */ \
    VAR(Float_t, met) /* type 1 corrected PF MET */  \
    VAR(Float_t, metphi) /* type 1 corrected PF MET phi */  \
    VAR(Float_t, puppimet) /* Puppi corrrected MET */  \
    VAR(Float_t, puppimetphi) /* Puppi corrected MET phi */  \
    VAR(Float_t, mvamet) /* mva met */  \
    VAR(Float_t, mvametphi) /* mva met phi */  \
    VAR(Float_t, pzetavis) /* see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pzetamiss) /* use MVA met, see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pfpzetamiss) /* As above but using pf met */  \
    VAR(Float_t, puppipzetamiss) /* As above but using puppi met */  \
    VAR(Float_t, mvacov00) /* mva met */  \
    VAR(Float_t, mvacov01) /* mva met */  \
    VAR(Float_t, mvacov10) /* mva met */  \
    VAR(Float_t, mvacov11) /* mva met */  \
    VAR(Float_t, metcov00) /* pf met */  \
    VAR(Float_t, metcov01) /* pf met */  \
    VAR(Float_t, metcov10) /* pf met */  \
    VAR(Float_t, metcov11) /* pf met */  \
    /* VBF system (Only fill if njetspt20>=2) */ \
    VAR(Float_t, mjj) /* (jet_1->vector()+jet_2->vector() ).M() */ \
    VAR(Float_t, jdeta) /* delta eta between leading jet and subleading jet */ \
    VAR(Float_t, njetingap) /* Number of jets passing pfJetID and pt > 30 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, njetingap20) /* Number of jets passing pfJetID and pt > 20 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, jdphi) /* delta phi between leading jet and subleading jet */ \
    /* additional jets */ \
    VAR(Int_t, nbtag) /* pt>20 and abs(eta)<2.4 */ \
    VAR(Int_t, njets) /* pt>30 and abs(eta)<4.7 */ \
    VAR(Int_t, njetspt20) /* pt>20 and abs(eta)<4.7 */ \
    JET_DATA(j, 1) /* leading jet sorted by pt (Fill only if corrected jet pt > 20 GeV) */ \
    JET_DATA(j, 2) /* trailing jet sorted by pt (Fill only if corrected jet pt>20 GeV) */ \
    JET_DATA(b, 1) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(b, 2) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    /* Extra lepton vetos */ \
    VAR(Bool_t, dilepton_veto) /* Event is vetoed by the dilepton veto if true */ \
    VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    VAR(Float_t, puweight) \
    /* hh->bbtautau part */ \
    VAR(Int_t, nbjets) /* pt>30 and abs(eta)<2.4 */ \
    JET_DATA(bjet_, 1) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(bjet_, 2) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    VAR(Float_t, m_kinfit) \
    VAR(Int_t, kinfit_convergence) \
    VAR(Float_t, deltaR_ll) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(htt_sync, SyncEvent, SyncTuple, SYNC_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(htt_sync, SyncTuple, SYNC_DATA)
#undef VAR
#undef SYNC_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef JVAR
