/*! Definiton of ntuple::SyncTree and ntuple::Sync classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define SYNC_DATA() \
    SIMPLE_VAR(Int_t, run) /* Run */ \
    SIMPLE_VAR(Int_t, lumi) /* Lumi */ \
    SIMPLE_VAR(Int_t, evt) /* Evt */ \
    /* Event Variables */ \
    SIMPLE_VAR(Int_t, npv) /* NPV */ \
    SIMPLE_VAR(Int_t, npu) /* NPU */ \
    SIMPLE_VAR(Float_t, rho) /* Rho */ \
    /* Event Weights */ \
    SIMPLE_VAR(Float_t, mcweight) /* MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs)) */ \
    SIMPLE_VAR(Float_t, puweight) /* Pielup Weight */ \
    SIMPLE_VAR(Float_t, trigweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, trigweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, idweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, idweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, isoweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, isoweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, fakeweight) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, effweight) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, DYweight) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, decayModeWeight_1) /* Tau decay mode weight for the first leg. */ \
    SIMPLE_VAR(Float_t, decayModeWeight_2) /* Tau decay mode weight for the second leg. */ \
    SIMPLE_VAR(Float_t, weight) /* mcweight*puweight*effweight */ \
    SIMPLE_VAR(Float_t, embeddedWeight) /*  */ \
    SIMPLE_VAR(Float_t, signalWeight) /*  */ \
    SIMPLE_VAR(Float_t, etau_fakerate) /* e to tau fake weight*/ \
    /* SV Fit variables */ \
    SIMPLE_VAR(Float_t, mvis) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, m_sv) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, pt_sv) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, eta_sv) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, phi_sv) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, m_sv_Up) /* High Energy scale shape */ \
    SIMPLE_VAR(Float_t, m_sv_Down) /* Low Energy Scale Shape */ \
    /* First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau */ \
    SIMPLE_VAR(Float_t, pt_1) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1) /* Eta */ \
    SIMPLE_VAR(Float_t, m_1) /* Mass */ \
    SIMPLE_VAR(Int_t, q_1) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_1) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_1) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_1) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, iso_1) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    SIMPLE_VAR(Float_t, mva_1) /* MVA id (when using electron) 0 otherwise */ \
    SIMPLE_VAR(Bool_t, passid_1) /* Whether it passes id  (not necessarily iso) */ \
    SIMPLE_VAR(Bool_t, passiso_1) /* Whether it passes iso (not necessarily id) */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMVA3raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose2_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonMedium2_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight2_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    /* Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau */ \
    SIMPLE_VAR(Float_t, pt_2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2) /* Eta */ \
    SIMPLE_VAR(Float_t, m_2) /* Mass */ \
    SIMPLE_VAR(Int_t, q_2) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_2) /* mT of 2nd lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_2) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_2) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, iso_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, mva_2) /* MVA id (for anti electron id) */ \
    SIMPLE_VAR(Bool_t, passid_2) /* Whether it passes id  (not necessarily iso) */ \
    SIMPLE_VAR(Bool_t, passiso_2) /* Whether it passes iso (not necessarily id) */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMVA3raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose2_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonMedium2_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight2_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    /* Di-lepton */ \
    SIMPLE_VAR(Float_t, pt_tt) /* pT */ \
    /* Met related variables */ \
    SIMPLE_VAR(Float_t, met) /* pfmet */ \
    SIMPLE_VAR(Float_t, metphi) /* pfmet Phi */ \
    SIMPLE_VAR(Float_t, mvamet) /* mvamet */ \
    SIMPLE_VAR(Float_t, mvametphi) /* mvamet Phi */ \
    SIMPLE_VAR(Float_t, pzetavis) /* pZeta Visible */ \
    SIMPLE_VAR(Float_t, pzetamiss) /* pZeta Missing */ \
    /* MET covariance matrices */ \
    SIMPLE_VAR(Float_t, metcov00) /* pf met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, metcov01) /* pf met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, metcov10) /* pf met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, metcov11) /* pf met covariance matrix 11 */ \
    /* MVAMet covariance matrices */ \
    SIMPLE_VAR(Float_t, mvacov00) /* mva met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, mvacov01) /* mva met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, mvacov10) /* mva met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, mvacov11) /* mva met covariance matrix 11 */ \
    /* number of jets passing jet id ( pt > 30 ) */ \
    SIMPLE_VAR(Int_t, njets) /*  */ \
    SIMPLE_VAR(Int_t, njetspt20) /*  */ \
    /* First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau) */ \
    SIMPLE_VAR(Float_t, jpt_1) /* Jet Pt after corrections */ \
    SIMPLE_VAR(Float_t, jeta_1) /* Jet Eta */ \
    SIMPLE_VAR(Float_t, jphi_1) /* Jet Phi */ \
    SIMPLE_VAR(Float_t, jptraw_1) /* Jet Raw Pt (before corrections) */ \
    SIMPLE_VAR(Float_t, jptunc_1) /* Jet Unc (relative to Jet corrected pT) */ \
    SIMPLE_VAR(Float_t, jmva_1) /* Jet MVA id value */ \
    SIMPLE_VAR(Float_t, jlrm_1) /* Jet MVA id value */ \
    SIMPLE_VAR(Int_t, jctm_1) /* Jet MVA id value */ \
    SIMPLE_VAR(Bool_t, jpass_1) /* Whether Jet pass PU Id Loose WP */ \
    /* Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau) */ \
    SIMPLE_VAR(Float_t, jpt_2) /* Jet Pt after corrections */ \
    SIMPLE_VAR(Float_t, jeta_2) /* Jet Eta */ \
    SIMPLE_VAR(Float_t, jphi_2) /* Jet Phi */ \
    SIMPLE_VAR(Float_t, jptraw_2) /* Jet Raw Pt (before corrections) */ \
    SIMPLE_VAR(Float_t, jptunc_2) /* Jet Unc (relative to Jet corrected pT) */ \
    SIMPLE_VAR(Float_t, jmva_2) /* Jet MVA id value */ \
    SIMPLE_VAR(Float_t, jlrm_2) /* Jet MVA id value */ \
    SIMPLE_VAR(Int_t, jctm_2) /* Jet MVA id value */ \
    SIMPLE_VAR(Bool_t, jpass_2) /* Whether jet passes PU Id Loose WP */ \
    /* number of btags passing btag id (medium CSV WP) ( pt > 20 ) */ \
    SIMPLE_VAR(Int_t, nbtag) /*  */ \
    /* Candidate B Jets : leading jet (in CSV ordering) passing (pt > 20 + eta < 2.4) */ \
    SIMPLE_VAR(Float_t, bpt_1) /* Corrected BTag Pt */ \
    SIMPLE_VAR(Float_t, beta_1) /* Btag Eta */ \
    SIMPLE_VAR(Float_t, bphi_1) /* Btag Phi */ \
    SIMPLE_VAR(Float_t, bcsv_1) /* Btag CSV */ \
    /* Candidate B Jets : subleading jet (in CSV ordering) passing (pt > 20 + eta < 2.4) */ \
    SIMPLE_VAR(Float_t, bpt_2) /* Corrected BTag Pt */ \
    SIMPLE_VAR(Float_t, beta_2) /* Btag Eta */ \
    SIMPLE_VAR(Float_t, bphi_2) /* Btag Phi */ \
    SIMPLE_VAR(Float_t, bcsv_2) /* Btag CSV */ \
    /* Candidate B Jets : third jet (in CSV ordering) passing (pt > 20 + eta < 2.4) */ \
    SIMPLE_VAR(Float_t, bpt_3) /* Corrected BTag Pt */ \
    SIMPLE_VAR(Float_t, beta_3) /* Btag Eta */ \
    SIMPLE_VAR(Float_t, bphi_3) /* Btag Phi */ \
    SIMPLE_VAR(Float_t, bcsv_3) /* Btag CSV */ \
    /**/ \
    SIMPLE_VAR(Float_t, m_bb) /* Corrected BTag Pt */ \
    SIMPLE_VAR(Float_t, m_ttbb) /* Btag Eta */ \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Sync, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, SyncTree, SYNC_DATA, Sync, "sync", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, SyncTree, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef SYNC_DATA

namespace ntuple {
inline double DefaultFillValueForSyncTree() { return -10000; }
}
