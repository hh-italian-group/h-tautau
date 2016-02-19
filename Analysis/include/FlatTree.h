/*! Definiton of ntuple::FlatTree and ntuple::Flat classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include <limits>
#include "AnalysisTools/Core/include/SmartTree.h"

#define FLAT_DATA() \
    /* Event Variables */ \
    SIMPLE_VAR(Int_t, run) /* Run */ \
    SIMPLE_VAR(Int_t, lumi) /* Lumi */ \
    SIMPLE_VAR(Int_t, evt) /* Event */ \
    SIMPLE_VAR(Int_t, channel) /* Analysis channel as defined in analysis::Channel */ \
    SIMPLE_VAR(Int_t, eventEnergyScale) /* identifier of the applied energy scale */ \
    SIMPLE_VAR(Int_t, eventType) /* event type category */ \
    \
    \
    /* First signal lepton :  muon for MuTau, electron for ETau, leading (in pT) tau for TauTau */ \
    SIMPLE_VAR(Float_t, pt_1) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1) /* Eta */ \
    SIMPLE_VAR(Float_t, m_1) /* Mass */ \
    SIMPLE_VAR(Float_t, energy_1) /* Energy */ \
    SIMPLE_VAR(Int_t, q_1) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_1) /* mT of the first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_1) /* d0 with respect to the primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_1) /* dZ with respect to the primary vertex */ \
    /* Gen particle quantities of the first signal lepton matched with the truth */ \
    SIMPLE_VAR(Int_t, pdgId_1_MC) /* PDG ID or particles::NONEXISTENT, if there is no matched genParticle. */ \
    SIMPLE_VAR(Float_t, pt_1_MC) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1_MC) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1_MC) /* Eta */ \
    SIMPLE_VAR(Float_t, m_1_MC) /* Mass */ \
    SIMPLE_VAR(Float_t, pt_1_visible_MC) /* Visible pT */ \
    SIMPLE_VAR(Float_t, phi_1_visible_MC) /* Visible phi */ \
    SIMPLE_VAR(Float_t, eta_1_visible_MC) /* Visible eta */ \
    SIMPLE_VAR(Float_t, m_1_visible_MC) /* Visible mass */ \
    /* First lepton - electron & muon specific */ \
    SIMPLE_VAR(Float_t, pfRelIso_1) /* Delta Beta for muon and electron */ \
    SIMPLE_VAR(Float_t, mva_1) /* MVA id when using electron, 0 otherwise */ \
    SIMPLE_VAR(Bool_t, passid_1) /* Whether it passes id (not necessarily iso) */ \
    SIMPLE_VAR(Bool_t, passiso_1) /* Whether it passes iso (not necessarily id) */ \
    /* First lepton - hadronic tau specific */ \
    SIMPLE_VAR(Int_t, decayMode_1) /* tau decay mode as defined in ntuple::tau_id::hadronicDecayMode */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) /* tau raw isolation value */ \
    SIMPLE_VAR(Float_t, iso_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Bool_t, againstElectronLooseMVA_1) /* Whether tau passes loose MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMediumMVA_1) /* Whether tau passes medium MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTightMVA_1) /* Whether tau passes tight MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronVTightMVA_1) /* Whether tau passes very tight MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronLooseMVA_custom_1) /* Whether tau passes loose MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMediumMVA_custom_1) /* Whether tau passes medium MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTightMVA_custom_1) /* Whether tau passes tight MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronVTightMVA_custom_1) /* Whether tau passes very tight MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronLoose_1) /* Whether tau passes loose against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMedium_1) /* Whether tau passes medium against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTight_1) /* Whether tau passes tight against electron discriminator */ \
    SIMPLE_VAR(Float_t, againstElectronMVA3raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Bool_t, againstMuonLoose_1) /* Whether tau passes loose against muon discriminator */ \
    SIMPLE_VAR(Bool_t, againstMuonMedium_1) /* Whether tau passes medium against muon discriminator */ \
    SIMPLE_VAR(Bool_t, againstMuonTight_1) /* Whether tau passes tight against muon discriminator */ \
    \
    \
    /* Second lepton : hadronic tau for MuTau & ETau, trailing (in pT) tau for TauTau */ \
    SIMPLE_VAR(Float_t, pt_2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2) /* Eta */ \
    SIMPLE_VAR(Float_t, m_2) /* Mass */ \
    SIMPLE_VAR(Float_t, energy_2) /* Energy */ \
    SIMPLE_VAR(Int_t, q_2) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_2) /* mT of second lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_2) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_2) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, iso_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    /* Gen particle quantities of second signal lepton matched with the truth */ \
    SIMPLE_VAR(Int_t, pdgId_2_MC) /* PDG ID or particles::NONEXISTENT, if there is no matched genParticle. */ \
    SIMPLE_VAR(Float_t, pt_2_MC) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2_MC) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2_MC) /* Eta */ \
    SIMPLE_VAR(Float_t, m_2_MC) /* Mass */ \
    SIMPLE_VAR(Float_t, pt_2_visible_MC) /* Visible pT */ \
    SIMPLE_VAR(Float_t, phi_2_visible_MC) /* Visible phi */ \
    SIMPLE_VAR(Float_t, eta_2_visible_MC) /* Visible eta */ \
    SIMPLE_VAR(Float_t, m_2_visible_MC) /* Visible mass */ \
    /* Second lepton - hadronic tau specific */ \
    SIMPLE_VAR(Int_t, decayMode_2) /* tau decay mode as defined in ntuple::tau_id::hadronicDecayMode */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) /* tau raw isolation value */ \
    SIMPLE_VAR(Bool_t, againstElectronLooseMVA_2) /* Whether tau passes loose MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMediumMVA_2) /* Whether tau passes medium MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTightMVA_2) /* Whether tau passes tight MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronVTightMVA_2) /* Whether tau passes very tight MVA against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronLooseMVA_custom_2) /* Whether tau passes loose MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMediumMVA_custom_2) /* Whether tau passes medium MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTightMVA_custom_2) /* Whether tau passes tight MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronVTightMVA_custom_2) /* Whether tau passes very tight MVA against electron custom discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronLoose_2) /* Whether tau passes loose against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronMedium_2) /* Whether tau passes medium against electron discriminator */ \
    SIMPLE_VAR(Bool_t, againstElectronTight_2) /* Whether tau passes tight against electron discriminator */ \
    SIMPLE_VAR(Float_t, againstElectronMVA3raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Bool_t, againstMuonLoose_2) /* Whether tau passes loose against muon discriminator */ \
    SIMPLE_VAR(Bool_t, againstMuonMedium_2) /* Whether tau passes medium against muon discriminator */ \
    SIMPLE_VAR(Bool_t, againstMuonTight_2) /* Whether tau passes tight against muon discriminator */ \
    \
    \
    /* H_tautau variables */ \
    SIMPLE_VAR(Float_t, DeltaR_leptons) /* DeltaR between two legs of H_tautau candidate */ \
    SIMPLE_VAR(Float_t, mvis) /* Visible mass of H_tautau */ \
    SIMPLE_VAR(Float_t, m_sv_vegas) /* Mass of H_tautau corrected by svFit using integration method VEGAS */ \
    SIMPLE_VAR(Float_t, m_sv_MC) /* Mass of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, pt_sv_MC) /* Pt of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, eta_sv_MC) /* Eta of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, phi_sv_MC) /* Phi of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, pt_tt) /* pt of two legs of H_tautau without MVAMET */ \
    SIMPLE_VAR(Float_t, pt_tt_MET) /* pt of two legs of H_tautau with MVAMET */ \
    \
    \
    /* Kinematic fit variables */ \
    SIMPLE_VAR(Float_t, kinfit_bb_tt_mass) /* Four body mass calculated using kinematic fit */ \
    SIMPLE_VAR(Int_t, kinfit_bb_tt_convergence) /* Convergence of four body kinematic fit */ \
    SIMPLE_VAR(Float_t, kinfit_bb_tt_chi2) /* Chi-square of four body kinematic fit */ \
    SIMPLE_VAR(Float_t, kinfit_bb_tt_pull_balance) /* Pull balance of four body kinematic fit */ \
    \
    \
    /* Met related variables */ \
    SIMPLE_VAR(Float_t, met) /* pfmet */ \
    SIMPLE_VAR(Float_t, metphi) /* pfmet Phi */ \
    SIMPLE_VAR(Float_t, mvamet) /* mvamet */ \
    SIMPLE_VAR(Float_t, mvametphi) /* mvamet Phi */ \
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
    \
    \
    /* Useful info at gen level */ \
    SIMPLE_VAR(Int_t, pdgId_resonance_MC) /* PDG ID of X/MSSM_H or particles::NONEXISTENT, if it not present */ \
    SIMPLE_VAR(Float_t, pt_resonance_MC) /* pt of X/MSSM_H */ \
    SIMPLE_VAR(Float_t, eta_resonance_MC) /* eta of X/MSSM_H */ \
    SIMPLE_VAR(Float_t, phi_resonance_MC) /* phi of X/MSSM_H */ \
    SIMPLE_VAR(Float_t, mass_resonance_MC) /* mass of X/MSSM_H */ \
    SIMPLE_VAR(Int_t, pdgId_Htt_MC) /* PDG ID of H_tt or particles::NONEXISTENT, if it not present */ \
    SIMPLE_VAR(Float_t, pt_Htt_MC) /* pt of Htt */ \
    SIMPLE_VAR(Float_t, eta_Htt_MC) /* eta of Htt */ \
    SIMPLE_VAR(Float_t, phi_Htt_MC) /* phi of Htt */ \
    SIMPLE_VAR(Float_t, mass_Htt_MC) /* mass of Htt */ \
    SIMPLE_VAR(Int_t, pdgId_Hbb_MC) /* PDG ID of H_bb or particles::NONEXISTENT, if it not present */ \
    SIMPLE_VAR(Float_t, pt_Hbb_MC) /* pt of Hbb */ \
    SIMPLE_VAR(Float_t, eta_Hbb_MC) /* eta of Hbb */ \
    SIMPLE_VAR(Float_t, phi_Hbb_MC) /* phi of Hbb */ \
    SIMPLE_VAR(Float_t, mass_Hbb_MC) /* mass of Hbb */ \
    SIMPLE_VAR(Int_t, n_extraJets_MC) /* number extra jets */ \
    \
    \
    /* Jets info */ \
    SIMPLE_VAR(Int_t, njets) /* number of jets passing jet id ( pt > 30 ) */ \
    SIMPLE_VAR(Int_t, njetspt20) /* number of jets passing jet id ( pt > 20 ) */ \
    /* All jets pt > 30 sorted in pt after applying Jet energy corrections (excluding hadronic Tau) */ \
    VECTOR_VAR(Float_t, pt_jets) /* Jets Pt after corrections */ \
    VECTOR_VAR(Float_t, eta_jets) /* Jets Eta */ \
    VECTOR_VAR(Float_t, phi_jets) /* Jets Phi */ \
    VECTOR_VAR(Float_t, ptraw_jets) /* Jets Raw Pt (before corrections) */ \
    VECTOR_VAR(Float_t, ptunc_jets) /* Jet Unc (relative to Jet corrected pT) */ \
    VECTOR_VAR(Float_t, mva_jets) /* Jet MVA id value */ \
    VECTOR_VAR(Bool_t, passPU_jets) /* Whether Jet pass PU Id Loose WP */ \
    /* b-jets info */ \
    SIMPLE_VAR(Int_t, nBjets) /* number of btags not passing btag id (medium CSV WP) ( pt > 20 ) without re-tag applied */ \
    SIMPLE_VAR(Int_t, nBjets_retagged) /* number of btags passing btag id (medium CSV WP) ( pt > 20 ) with re-tag applied */ \
    /* All b-jets passing jet id ( pt > 20 ) sorted by CSV without re-tag applied */ \
    VECTOR_VAR(Float_t, pt_Bjets) /* pt of b-jets */ \
    VECTOR_VAR(Float_t, eta_Bjets) /* eta of b-jets */ \
    VECTOR_VAR(Float_t, phi_Bjets) /* phi of b-jets */ \
    VECTOR_VAR(Float_t, energy_Bjets) /* energy of b-jets */ \
    VECTOR_VAR(Float_t, chargedHadronEF_Bjets) /* Charged hadron energy fraction of b-jets */ \
    VECTOR_VAR(Float_t, neutralHadronEF_Bjets) /* Neutral hadron energy fraction of b-jets */ \
    VECTOR_VAR(Float_t, photonEF_Bjets) /* photon energy fraction of b-jets */ \
    VECTOR_VAR(Float_t, muonEF_Bjets) /* muon energy fraction of b-jets */ \
    VECTOR_VAR(Float_t, electronEF_Bjets) /* electron energy fraction of b-jets */ \
    VECTOR_VAR(Float_t, csv_Bjets) /* csv of b-jets */ \
    VECTOR_VAR(Bool_t, isBjet_MC_Bjet) /* Whether b-jet matches MC b-jet */ \
    VECTOR_VAR(Bool_t, isBjet_MC_Bjet_withLeptonicDecay) /* Whether b-jet matches MC b-jet that decayed leptonically */ \
    \
    \
    /* vertices */ \
    SIMPLE_VAR(Float_t, x_PV) /* x of PV*/ \
    SIMPLE_VAR(Float_t, y_PV) /* y of PV*/ \
    SIMPLE_VAR(Float_t, z_PV) /* z of PV*/ \
    SIMPLE_VAR(Int_t, npv) /* NPV */ \
    SIMPLE_VAR(Int_t, npu) /* NPU */ \
    /* Event Weights */ \
    SIMPLE_VAR(Float_t, puweight) /* Pielup weight */ \
    SIMPLE_VAR(Float_t, trigweight_1) /* Trigger weight for the first leg */ \
    SIMPLE_VAR(Float_t, trigweight_2) /* Trigger weight for the second leg */ \
    SIMPLE_VAR(Float_t, idweight_1) /* ID weight for the first leg */ \
    SIMPLE_VAR(Float_t, idweight_2) /* ID weight for the second leg */ \
    SIMPLE_VAR(Float_t, isoweight_1) /* Isolation weight for the first leg */ \
    SIMPLE_VAR(Float_t, isoweight_2) /* Isolation weight for the second leg */ \
    SIMPLE_VAR(Float_t, fakeweight_1) /* fake rate weight for the first leg (only e->tau)*/ \
    SIMPLE_VAR(Float_t, fakeweight_2) /* fake rate weight for the second leg */ \
    SIMPLE_VAR(Float_t, decayModeWeight_1) /* decay mode weight for the first leg */ \
    SIMPLE_VAR(Float_t, decayModeWeight_2) /* decay mode weight for the second leg */ \
    SIMPLE_VAR(Float_t, embeddedWeight) /* Weight for embedded events */ \
    SIMPLE_VAR(Float_t, weight) /* Product of all weights defined above */ \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Flat, FLAT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, FlatTree, FLAT_DATA, Flat, "flat", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, FlatTree, FLAT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef FLAT_DATA

namespace ntuple {
inline float DefaultFloatFillValueForFlatTree() { return std::numeric_limits<float>::lowest(); }
inline int DefaultIntegerFillValueForFlatTree() { return std::numeric_limits<int>::lowest(); }

enum class EventType { Unknown = 0, ZL = 1, ZJ = 2, ZTT = 3, ZTT_L = 4 };

namespace detail {
std::map<EventType, std::string> EventTypeNameMap = {
    { EventType::Unknown, "Unknown" }, { EventType::ZL, "ZL" }, { EventType::ZJ, "ZJ" }, { EventType::ZTT, "ZTT" },
    { EventType::ZTT_L, "ZTT_L" }
};
} // namespace detail

std::ostream& operator<< (std::ostream& s, const EventType& t)
{
    s << detail::EventTypeNameMap.at(t);
    return s;
}

} // namespace ntuple
