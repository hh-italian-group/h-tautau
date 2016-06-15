/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define EVENT_DATA() \
    VAR(Int_t, run) /* Run */ \
    VAR(Int_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    VAR(Int_t, channelID) /* Channel: MuTau, ETau, TauTau */ \
    VAR(Int_t, eventEnergyScale) /* event type category */ \
    VAR(Int_t, eventType) /* event type category */ \
    VAR(Int_t, HTBin) /* event type category */ \
    VAR(Double_t, weightevt) /*Gen Event Weight*/ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* NPV */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(Float_t, rho) /* Use fixedGridRhoFastjetAll */ \
    /* SV Fit variables */ \
    VAR(Float_t, m_sv) /* SV Fit using integration method */ \
    VAR(Float_t, pt_sv) /* SV Fit using integration method */ \
    VAR(Float_t, eta_sv) /* SV Fit using integration method */ \
    VAR(Float_t, phi_sv) /* SV Fit using integration method */ \
    /* First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau */ \
    VAR(Float_t, pt_1) /* pT */ \
    VAR(Float_t, phi_1) /* Phi */ \
    VAR(Float_t, eta_1) /* Eta */ \
    VAR(Float_t, m_1) /* Mass */ \
    VAR(Int_t, q_1) /* Charge */ \
    VAR(Float_t, d0_1) /* d0 with respect to primary vertex */ \
    VAR(Float_t, dZ_1) /* dZ with respect to primary vertex */ \
    VAR(Float_t, mt_1) /* mT of  first lepton wrt to MVA met */ \
    VAR(Float_t, pfmt_1) /* mT of  first lepton wrt to MVA met */ \
    VAR(Float_t, puppimt_1) /* mT of  first lepton wrt to MVA met */ \
    VAR(Float_t, iso_1) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    VAR(Float_t, id_e_mva_nt_loose_1) /* Non-triggering electron ID MVA score id (when using electron) 0 otherwise */ \
    VAR(Int_t, gen_match_1 ) /*Generator matching, see Htautau Twiki*/\
    VAR(Float_t, againstElectronLooseMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronMediumMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronTightMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronVLooseMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronVTightMVA6_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, byIsolationMVA2raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstMuonLoose3_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstMuonTight3_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) /*  */ \
    VAR(Float_t, byIsolationMVA3newDMwoLTraw_1) /* MVA iso for Tau w/o Lifetime information New Decay Mode */ \
    VAR(Float_t, byIsolationMVA3oldDMwoLTraw_1) /* MVA iso for Tau w/o Lifetime information Old Decay Mode */ \
    VAR(Float_t, byIsolationMVA3newDMwLTraw_1) /* MVA iso for Tau w/ Lifetime information New Decay Mode */ \
    VAR(Float_t, byIsolationMVA3oldDMwLTraw_1) /* MVA iso for Tau w/ Lifetime information Old Decay Mode */ \
    VAR(Int_t, byVLooseIsolationMVArun2v1DBoldDMwLT_1) /* MVA iso for Tau w/ Lifetime VLoose WP */ \
    VAR(Int_t, byLooseIsolationMVArun2v1DBoldDMwLT_1) /* MVA iso for Tau w/ Lifetime Loose WP */ \
    VAR(Int_t, byMediumIsolationMVArun2v1DBoldDMwLT_1) /* MVA iso for Tau w/ Lifetime Medium WP */ \
    VAR(Int_t, byTightIsolationMVArun2v1DBoldDMwLT_1) /* MVA iso for Tau w/ Lifetime Tight WP */ \
    VAR(Int_t, byVTightIsolationMVArun2v1DBoldDMwLT_1) /* MVA iso for Tau w/ Lifetime VTight WP */ \
    VAR(Float_t, chargedIsoPtSum_1) \
    VAR(Int_t, decayModeFindingOldDMs_1) /* Old Decay Mode finding */\
    VAR(Float_t, neutralIsoPtSum_1) \
    VAR(Float_t, puCorrPtSum_1) \
    VAR(Float_t, trigweight_1) \
    VAR(Float_t, idisoweight_1) \
    /* Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau */ \
    VAR(Float_t, pt_2) /* pT */ \
    VAR(Float_t, phi_2) /* Phi */ \
    VAR(Float_t, eta_2) /* Eta */ \
    VAR(Float_t, m_2) /* Mass */ \
    VAR(Int_t, q_2) /* Charge */ \
    VAR(Float_t, d0_2) /* d0 with respect to primary vertex */ \
    VAR(Float_t, dZ_2) /* dZ with respect to primary vertex */ \
    VAR(Float_t, mt_2) /* mT of  first lepton wrt to PF met */ \
    VAR(Float_t, pfmt_2) /* mT of  first lepton wrt to PF met */ \
    VAR(Float_t, puppimt_2) /* mT of  first lepton wrt to PF met */ \
    VAR(Float_t, iso_2) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    VAR(Int_t, gen_match_2 ) /*Generator matching, see Htautau Twiki*/\
    VAR(Float_t, againstElectronLooseMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronMediumMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronTightMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronVLooseMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstElectronVTightMVA6_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, byIsolationMVA2raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstMuonLoose3_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, againstMuonTight3_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) /*  */ \
    VAR(Float_t, byIsolationMVA3newDMwoLTraw_2) /* MVA iso for Tau w/o Lifetime information New Decay Mode */ \
    VAR(Float_t, byIsolationMVA3oldDMwoLTraw_2) /* MVA iso for Tau w/o Lifetime information Old Decay Mode */ \
    VAR(Float_t, byIsolationMVA3newDMwLTraw_2) /* MVA iso for Tau w/ Lifetime information New Decay Mode */ \
    VAR(Float_t, byIsolationMVA3oldDMwLTraw_2) /* MVA iso for Tau w/ Lifetime information Old Decay Mode */ \
    VAR(Int_t, byVLooseIsolationMVArun2v1DBoldDMwLT_2) /* MVA iso for Tau w/ Lifetime VLoose WP */ \
    VAR(Int_t, byLooseIsolationMVArun2v1DBoldDMwLT_2) /* MVA iso for Tau w/ Lifetime Loose WP */ \
    VAR(Int_t, byMediumIsolationMVArun2v1DBoldDMwLT_2) /* MVA iso for Tau w/ Lifetime Medium WP */ \
    VAR(Int_t, byTightIsolationMVArun2v1DBoldDMwLT_2) /* MVA iso for Tau w/ Lifetime Tight WP */ \
    VAR(Int_t, byVTightIsolationMVArun2v1DBoldDMwLT_2) /* MVA iso for Tau w/ Lifetime VTight WP */ \
    VAR(Float_t, chargedIsoPtSum_2) \
    VAR(Int_t, decayModeFindingOldDMs_2) /* Old Decay Mode finding */\
    VAR(Float_t, neutralIsoPtSum_2) \
    VAR(Float_t, puCorrPtSum_2) \
    VAR(Float_t, idisoweight_2) \
    /* Di-lepton */ \
    VAR(Float_t, m_vis) /* pairs invariant mass */ \
    VAR(Float_t, pt_tt) /* pT */ \
    /* Met related variables */ \
    VAR(Float_t, met) /* pfmet */ \
    VAR(Float_t, metphi) /* pfmet Phi */ \
    VAR(Float_t, puppimet) /* puppimet */ \
    VAR(Float_t, puppimetphi) /* puppimet Phi */ \
    VAR(Bool_t, isPFMET) /* pfmet Phi */ \
    VAR(Float_t, mvamet) /* mvamet */ \
    VAR(Float_t, mvametphi) /* mvamet Phi */ \
    VAR(Float_t, pzetavis) /* pZeta Visible */ \
    VAR(Float_t, pzetamiss) /* pZeta Missing */ \
    /* MVAMet covariance matrices */ \
    VAR(Float_t, mvacov00) /* mva met covariance matrix 00 */ \
    VAR(Float_t, mvacov01) /* mva met covariance matrix 01 */ \
    VAR(Float_t, mvacov10) /* mva met covariance matrix 10 */ \
    VAR(Float_t, mvacov11) /* mva met covariance matrix 11 */ \
    /* MET covariance matrices */ \
    VAR(Float_t, metcov00) /* pf met covariance matrix 00 */ \
    VAR(Float_t, metcov01) /* pf met covariance matrix 01 */ \
    VAR(Float_t, metcov10) /* pf met covariance matrix 10 */ \
    VAR(Float_t, metcov11) /* pf met covariance matrix 11 */ \
    /* number of jets passing jet id ( pt > 30 ) */ \
    VAR(Int_t, njets) /*  */ \
    VAR(Int_t, njetspt20) /*  */ \
    /* Candidate Jets: jets after applying Jet energy corrections (excluding hadronic Tau) */ \
    VAR(std::vector<Float_t>, pt_jets) /* Jet Pt after corrections */ \
    VAR(std::vector<Float_t>, eta_jets) /* Jet Eta */ \
    VAR(std::vector<Float_t>, phi_jets) /* Jet Phi */ \
    VAR(std::vector<Float_t>, energy_jets) /* Jet Energy */ \
    VAR(std::vector<Float_t>, rawf_jets) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    VAR(std::vector<Float_t>, mva_jets) /* Jet MVA id value */ \
    VAR(std::vector<Float_t>, csv_jets) /* Jet CSV value */ \
    VAR(std::vector<Int_t>, partonFlavour_jets) \
    /* number of btags passing btag id (medium CSV WP) ( pt > 20 ) */ \
    VAR(Int_t, nbtag) /*  */ \
    /* Candidate B Jets (in pt ordering) passing (pt > 20 + eta < 2.4) */ \
    VAR(std::vector<Float_t>, pt_bjets) /* Corrected BTag Pt */ \
    VAR(std::vector<Float_t>, eta_bjets) /* Btag Eta */ \
    VAR(std::vector<Float_t>, phi_bjets) /* Btag Phi */ \
    VAR(std::vector<Float_t>, energy_bjets) /* Btag Energy */ \
    VAR(std::vector<Float_t>, rawf_bjets) /* Btag factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    VAR(std::vector<Float_t>, mva_bjets) /* Btag mva */ \
    VAR(std::vector<Float_t>, csv_bjets) /* Btag CSV */ \
    VAR(std::vector<Int_t>, partonFlavour_bjets) /* Jet CSV value */ \
    /* KinFit Variables */ \
    VAR(std::vector<Double_t>, kinFit_m) /* KinFit m_bbtt mass compute the first 2 jets, ordered by CSV*/\
    VAR(std::vector<Double_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Double_t>, kinFit_probability) /*  KinFit chi2 probability value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    /* LHE info */\
    VAR(std::vector<Int_t>, lhe_particle_pdg) \
    VAR(std::vector<Float_t>, lhe_particle_pt) \
    VAR(std::vector<Float_t>, lhe_particle_eta) \
    VAR(std::vector<Float_t>, lhe_particle_phi) \
    VAR(std::vector<Float_t>, lhe_particle_m) \
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
