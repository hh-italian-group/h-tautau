
#pragma once

//#include "TreeProduction/interface/SmartTree.h"
#include "../../TreeProduction/interface/SmartTree.h"

#define SYNC_DATA() \
    SIMPLE_VAR(Int_t, run) /* Run */ \
    SIMPLE_VAR(Int_t, lumi) /* Lumi */ \
    SIMPLE_VAR(ULong64_t, evt) /* Evt */ \
    SIMPLE_VAR(Int_t, eventType) /* event type category */ \
    SIMPLE_VAR(Int_t, HTBin) /* event type category */ \
    SIMPLE_VAR(Double_t, weightevt) /*Gen Event Weight*/ \
    /* Event Variables */ \
    SIMPLE_VAR(Int_t, npv) /* NPV */ \
    SIMPLE_VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    SIMPLE_VAR(Float_t, rho) /* Use fixedGridRhoFastjetAll */ \
    /* SV Fit variables */ \
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
    SIMPLE_VAR(Float_t, d0_1) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_1) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, mt_1) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, pfmt_1) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, puppimt_1) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, iso_1) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    SIMPLE_VAR(Float_t, id_e_mva_nt_loose_1) /* Non-triggering electron ID MVA score id (when using electron) 0 otherwise */ \
    SIMPLE_VAR(Int_t, gen_match_1 ) /*Generator matching, see Htautau Twiki*/\
    SIMPLE_VAR(Float_t, againstElectronLooseMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronMediumMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronTightMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronVLooseMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronVTightMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose3_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight3_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) /*  */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwoLTraw_1) /* MVA iso for Tau w/o Lifetime information New Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwoLTraw_1) /* MVA iso for Tau w/o Lifetime information Old Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwLTraw_1) /* MVA iso for Tau w/ Lifetime information New Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwLTraw_1) /* MVA iso for Tau w/ Lifetime information Old Decay Mode */ \
    SIMPLE_VAR(Float_t, chargedIsoPtSum_1) \
    SIMPLE_VAR(Int_t, decayModeFindingOldDMs_1) /* Old Decay Mode finding */\
    SIMPLE_VAR(Float_t, neutralIsoPtSum_1) \
    SIMPLE_VAR(Float_t, puCorrPtSum_1) \
    SIMPLE_VAR(Float_t, trigweight_1) \
    SIMPLE_VAR(Float_t, idisoweight_1) \
    /* Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau */ \
    SIMPLE_VAR(Float_t, pt_2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2) /* Eta */ \
    SIMPLE_VAR(Float_t, m_2) /* Mass */ \
    SIMPLE_VAR(Int_t, q_2) /* Charge */ \
    SIMPLE_VAR(Float_t, d0_2) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_2) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, mt_2) /* mT of  first lepton wrt to PF met */ \
    SIMPLE_VAR(Float_t, pfmt_2) /* mT of  first lepton wrt to PF met */ \
    SIMPLE_VAR(Float_t, puppimt_2) /* mT of  first lepton wrt to PF met */ \
    SIMPLE_VAR(Float_t, iso_2) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    SIMPLE_VAR(Int_t, gen_match_2 ) /*Generator matching, see Htautau Twiki*/\
    SIMPLE_VAR(Float_t, againstElectronLooseMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronMediumMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronTightMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronVLooseMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronVTightMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose3_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight3_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) /*  */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwoLTraw_2) /* MVA iso for Tau w/o Lifetime information New Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwoLTraw_2) /* MVA iso for Tau w/o Lifetime information Old Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwLTraw_2) /* MVA iso for Tau w/ Lifetime information New Decay Mode */ \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwLTraw_2) /* MVA iso for Tau w/ Lifetime information Old Decay Mode */ \
    SIMPLE_VAR(Float_t, chargedIsoPtSum_2) \
    SIMPLE_VAR(Int_t, decayModeFindingOldDMs_2) /* Old Decay Mode finding */\
    SIMPLE_VAR(Float_t, neutralIsoPtSum_2) \
    SIMPLE_VAR(Float_t, puCorrPtSum_2) \
    SIMPLE_VAR(Float_t, idisoweight_2) \
    /* Di-lepton */ \
    SIMPLE_VAR(Float_t, m_vis) /* pairs invariant mass */ \
    SIMPLE_VAR(Float_t, pt_tt) /* pT */ \
    /* Met related variables */ \
    SIMPLE_VAR(Float_t, met) /* pfmet */ \
    SIMPLE_VAR(Float_t, metphi) /* pfmet Phi */ \
    SIMPLE_VAR(Float_t, puppimet) /* puppimet */ \
    SIMPLE_VAR(Float_t, puppimetphi) /* puppimet Phi */ \
    SIMPLE_VAR(Bool_t, isPFMET) /* pfmet Phi */ \
    SIMPLE_VAR(Float_t, mvamet) /* mvamet */ \
    SIMPLE_VAR(Float_t, mvametphi) /* mvamet Phi */ \
    SIMPLE_VAR(Float_t, pzetavis) /* pZeta Visible */ \
    SIMPLE_VAR(Float_t, pzetamiss) /* pZeta Missing */ \
    /* MVAMet covariance matrices */ \
    SIMPLE_VAR(Float_t, mvacov00) /* mva met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, mvacov01) /* mva met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, mvacov10) /* mva met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, mvacov11) /* mva met covariance matrix 11 */ \
    /* MET covariance matrices */ \
    SIMPLE_VAR(Float_t, metcov00) /* pf met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, metcov01) /* pf met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, metcov10) /* pf met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, metcov11) /* pf met covariance matrix 11 */ \
    /* number of jets passing jet id ( pt > 30 ) */ \
    SIMPLE_VAR(Int_t, njets) /*  */ \
    SIMPLE_VAR(Int_t, njetspt20) /*  */ \
    /* Candidate Jets: jets after applying Jet energy corrections (excluding hadronic Tau) */ \
    VECTOR_VAR(Float_t, pt_jets) /* Jet Pt after corrections */ \
    VECTOR_VAR(Float_t, eta_jets) /* Jet Eta */ \
    VECTOR_VAR(Float_t, phi_jets) /* Jet Phi */ \
    VECTOR_VAR(Float_t, energy_jets) /* Jet Energy */ \
    VECTOR_VAR(Float_t, rawf_jets) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    VECTOR_VAR(Float_t, mva_jets) /* Jet MVA id value */ \
    VECTOR_VAR(Float_t, csv_jets) /* Jet CSV value */ \
    /* Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau) */ \
    /* SIMPLE_VAR(Float_t, jpt_2)  Jet Pt after corrections */ \
    /* SIMPLE_VAR(Float_t, jeta_2)  Jet Eta */ \
    /* SIMPLE_VAR(Float_t, jphi_2)  Jet Phi */ \
    /* SIMPLE_VAR(Float_t, jrawf_2)  factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    /* SIMPLE_VAR(Float_t, jmva_2)  Jet MVA id value */ \
    /* number of btags passing btag id (medium CSV WP) ( pt > 20 ) */ \
    SIMPLE_VAR(Int_t, nbtag) /*  */ \
    /* Candidate B Jets (in pt ordering) passing (pt > 20 + eta < 2.4) */ \
    VECTOR_VAR(Float_t, pt_bjets) /* Corrected BTag Pt */ \
    VECTOR_VAR(Float_t, eta_bjets) /* Btag Eta */ \
    VECTOR_VAR(Float_t, phi_bjets) /* Btag Phi */ \
    VECTOR_VAR(Float_t, energy_bjets) /* Btag Energy */ \
    VECTOR_VAR(Float_t, rawf_bjets) /* Btag factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    VECTOR_VAR(Float_t, mva_bjets) /* Btag mva */ \
    VECTOR_VAR(Float_t, csv_bjets) /* Btag CSV */ \
    /* Candidate B Jets : subleading jet (in CSV ordering) passing (pt > 20 + eta < 2.4) */ \
    /* SIMPLE_VAR(Float_t, bpt_2)  Corrected BTag Pt */ \
    /* SIMPLE_VAR(Float_t, beta_2)  Btag Eta */ \
    /* SIMPLE_VAR(Float_t, bphi_2)  Btag Phi */ \
    /* SIMPLE_VAR(Float_t, brawf_2)  Btag factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    /* SIMPLE_VAR(Float_t, bmva_2)  Btag mva */ \
    /* SIMPLE_VAR(Float_t, bcsv_2)  Btag CSV */ \
    /**/ \
    SIMPLE_VAR(Float_t, HT) \
    SIMPLE_VAR(Bool_t, dilepton_veto) /* Event is vetoed by the dilepton veto if true */ \
    SIMPLE_VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    SIMPLE_VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(Run2, Sync, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(Run2, SyncTree, SYNC_DATA, Sync, "sync", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(Run2, SyncTree, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef SYNC_DATA

namespace Run2 {
inline double DefaultFillValueForSyncTree() { return -10000; }
inline float DefaultFloatFillValueForSyncTree() { return std::numeric_limits<float>::lowest(); }

enum class HTbinning { lt100 = 0, f100to200 = 1, f200to400 = 2, f400to600 = 3, gt600 = 4 };

}


