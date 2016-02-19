/*! Definiton of ntuple::TauTree and ntuple::Tau classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define TAU_DATA() \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, mass) \
    /* Charge */ \
    SIMPLE_VAR(Int_t, charge) \
    /* Decay mode */ \
    SIMPLE_VAR(Int_t, decayMode) \
    /* Leading particle pT */ \
    SIMPLE_VAR(Float_t, leadChargedParticlePt) \
    SIMPLE_VAR(Float_t, leadNeutralParticleEt) \
    SIMPLE_VAR(Float_t, leadParticleEt) \
    /* Number of charged/neutral candidates and photons in different cones */ \
    SIMPLE_VAR(UInt_t, numChargedHadronsSignalCone) \
    SIMPLE_VAR(UInt_t, numNeutralHadronsSignalCone) \
    SIMPLE_VAR(UInt_t, numPhotonsSignalCone) \
    SIMPLE_VAR(UInt_t, numParticlesSignalCone) \
    SIMPLE_VAR(UInt_t, numChargedHadronsIsoCone) \
    SIMPLE_VAR(UInt_t, numNeutralHadronsIsoCone) \
    SIMPLE_VAR(UInt_t, numPhotonsIsoCone) \
    SIMPLE_VAR(UInt_t, numParticlesIsoCone) \
    SIMPLE_VAR(Float_t, ptSumPFChargedHadronsIsoCone) \
    SIMPLE_VAR(Float_t, etSumPhotonsIsoCone) \
    /* Charged hadron candidates */ \
    VECTOR_VAR(Float_t, signalChHadCand_Pt) \
    VECTOR_VAR(Float_t, signalChHadCand_Eta) \
    VECTOR_VAR(Float_t, signalChHadCand_Phi) \
    VECTOR_VAR(Float_t, isoChHadCand_Pt) \
    VECTOR_VAR(Float_t, isoChHadCand_Eta) \
    VECTOR_VAR(Float_t, isoChHadCand_Phi) \
    /* Neutral hadron candidates */ \
    VECTOR_VAR(Float_t, signalNeutrHadCand_Pt) \
    VECTOR_VAR(Float_t, signalNeutrHadCand_Eta) \
    VECTOR_VAR(Float_t, signalNeutrHadCand_Phi) \
    VECTOR_VAR(Float_t, isoNeutrHadCand_Pt) \
    VECTOR_VAR(Float_t, isoNeutrHadCand_Eta) \
    VECTOR_VAR(Float_t, isoNeutrHadCand_Phi) \
    /* Gamma candidates */ \
    VECTOR_VAR(Float_t, signalGammaCand_Pt) \
    VECTOR_VAR(Float_t, signalGammaCand_Eta) \
    VECTOR_VAR(Float_t, signalGammaCand_Phi) \
    VECTOR_VAR(Float_t, isoGammaCand_Pt) \
    VECTOR_VAR(Float_t, isoGammaCand_Eta) \
    VECTOR_VAR(Float_t, isoGammaCand_Phi) \
   /* MVA */ \
    SIMPLE_VAR(Float_t, leadPFCand_mva_e_pi) \
    /* PFTau specific variables DataFormats/PatCandidates/interface/TauPFSpecific.h */ \
    /* see DataFormats/TauReco/interface/PFTau.h */ \
    SIMPLE_VAR(Float_t, emFraction) \
    SIMPLE_VAR(Float_t, maximumHCALPFClusterEt) \
    SIMPLE_VAR(Float_t, ecalStripSumEOverPLead) \
    SIMPLE_VAR(Float_t, bremsRecoveryEOverPLead) \
    SIMPLE_VAR(Float_t, hcalTotOverPLead) \
    SIMPLE_VAR(Float_t, hcalMaxOverPLead) \
    SIMPLE_VAR(Float_t, hcal3x3OverPLead) \
    /* ET weighted eta-phi statistics starting from PFJet, see DataFormats/JetReco/interface/Jet.h */ \
    SIMPLE_VAR(Float_t, etaetaMoment) \
    SIMPLE_VAR(Float_t, phiphiMoment) \
    SIMPLE_VAR(Float_t, etaphiMoment)  \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vx) \
    SIMPLE_VAR(Float_t, vy) \
    SIMPLE_VAR(Float_t, vz) \
    /* Trigger match information */ \
    VECTOR_VAR(std::string, matchedTriggerPaths) \
    /* tau discriminators */ \
    TAU_DISCRIMINATOR_DATA() \
    /**/

// new TauID
#define NEW_TAU_DISCRIMINATOR_DATA() \
    /* Discriminator by decay mode */ \
    SIMPLE_VAR(Float_t, decayModeFinding) \
    /* Discriminators agains electron */ \
    SIMPLE_VAR(Float_t, againstElectronLoose) \
    SIMPLE_VAR(Float_t, againstElectronMedium) \
    SIMPLE_VAR(Float_t, againstElectronTight) \
    SIMPLE_VAR(Float_t, againstElectronLooseMVA5) \
    SIMPLE_VAR(Float_t, againstElectronMediumMVA5) \
    SIMPLE_VAR(Float_t, againstElectronTightMVA5) \
    SIMPLE_VAR(Float_t, againstElectronVTightMVA5) \
    /* Discriminators agains muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose) \
    SIMPLE_VAR(Float_t, againstMuonMedium) \
    SIMPLE_VAR(Float_t, againstMuonTight) \
    SIMPLE_VAR(Float_t, againstMuonLoose3) \
    SIMPLE_VAR(Float_t, againstMuonTight3) \
    SIMPLE_VAR(Float_t, againstMuonLooseMVA) \
    SIMPLE_VAR(Float_t, againstMuonMediumMVA) \
    SIMPLE_VAR(Float_t, againstMuonTightMVA) \
    SIMPLE_VAR(Float_t, againstMuonMVAraw) \
    /* Discriminators for isolation */ \
    SIMPLE_VAR(Float_t, byVLooseCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr3Hits) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr3Hits) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr3Hits) \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwoLTraw) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3oldDMwoLT) \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwLTraw) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3oldDMwLT) \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwoLTraw) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3newDMwoLT) \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwLTraw) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3newDMwLT) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3newDMwLT) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3newDMwLT) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3newDMwLT) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3newDMwLT) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3newDMwLT) \
//    /**/

// old TauID
#define TAU_DISCRIMINATOR_DATA() \
    /* Discriminator by decay mode */ \
    SIMPLE_VAR(Float_t, decayModeFinding) \
    /* Discriminators agains electron */ \
    SIMPLE_VAR(Float_t, againstElectronLoose) \
    SIMPLE_VAR(Float_t, againstElectronMedium) \
    SIMPLE_VAR(Float_t, againstElectronTight) \
    SIMPLE_VAR(Float_t, againstElectronMVA3raw) \
    SIMPLE_VAR(Float_t, againstElectronMVA3category) \
    SIMPLE_VAR(Float_t, againstElectronLooseMVA3) \
    SIMPLE_VAR(Float_t, againstElectronMediumMVA3) \
    SIMPLE_VAR(Float_t, againstElectronTightMVA3) \
    SIMPLE_VAR(Float_t, againstElectronVTightMVA3) \
    SIMPLE_VAR(Float_t, againstElectronDeadECAL) \
    /* Discriminators agains muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose) \
    SIMPLE_VAR(Float_t, againstMuonMedium) \
    SIMPLE_VAR(Float_t, againstMuonTight) \
    SIMPLE_VAR(Float_t, againstMuonLoose2) \
    SIMPLE_VAR(Float_t, againstMuonMedium2) \
    SIMPLE_VAR(Float_t, againstMuonTight2) \
    /* Discriminators for isolation */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw) \
    SIMPLE_VAR(Float_t, byVLooseCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr) \
    SIMPLE_VAR(Float_t, byIsolationMVAraw) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA) \
    SIMPLE_VAR(Float_t, byIsolationMVA2raw) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA2) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA2) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA2) \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr3Hits) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr3Hits) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr3Hits) \
    /**/


#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Tau, TAU_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, TauTree, TAU_DATA, Tau, "taus", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, TauTree, TAU_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TAU_DATA

namespace ntuple {
namespace tau_id {
enum hadronicDecayMode {
  kNull = -1,
  kOneProng0PiZero,
  kOneProng1PiZero,
  kOneProng2PiZero,
  kOneProng3PiZero,
  kOneProngNPiZero,
  kTwoProng0PiZero,
  kTwoProng1PiZero,
  kTwoProng2PiZero,
  kTwoProng3PiZero,
  kTwoProngNPiZero,
  kThreeProng0PiZero,
  kThreeProng1PiZero,
  kThreeProng2PiZero,
  kThreeProng3PiZero,
  kThreeProngNPiZero,
  kRareDecayMode
};

template<typename Value>
inline hadronicDecayMode ConvertToHadronicDecayMode(const Value& value)
{
    if(value < kNull || value > kRareDecayMode)
        throw std::runtime_error("value is not a hadronicDecayMode");
    return static_cast<ntuple::tau_id::hadronicDecayMode>(value);
}

}
}
