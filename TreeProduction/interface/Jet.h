/*! Definiton of ntuple::JetTree and ntuple::Jet classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define JET_DATA() \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, mass) \
    SIMPLE_VAR(Float_t, pt_raw) \
    SIMPLE_VAR(Float_t, energy_raw) \
    /* Jet energy correction */ \
    SIMPLE_VAR(Int_t, partonFlavour) \
    /* Jet identification in high PU environment */ \
    SIMPLE_VAR(Float_t, puIdMVA) \
    SIMPLE_VAR(Int_t, puIdBits) \
    SIMPLE_VAR(Float_t, puIdMVA_met) \
    SIMPLE_VAR(Int_t, puIdBits_met) \
    SIMPLE_VAR(Float_t, correction) \
    /* Energy fraction information */ \
    SIMPLE_VAR(Float_t, chargedEmEnergyFraction) \
    SIMPLE_VAR(Float_t, chargedHadronEnergyFraction) \
    SIMPLE_VAR(Float_t, chargedMuEnergyFraction) \
    SIMPLE_VAR(Float_t, electronEnergyFraction) \
    SIMPLE_VAR(Float_t, muonEnergyFraction) \
    SIMPLE_VAR(Float_t, neutralEmEnergyFraction) \
    SIMPLE_VAR(Float_t, neutralHadronEnergyFraction) \
    SIMPLE_VAR(Float_t, photonEnergyFraction) \
    /* Multiplicity information */ \
    SIMPLE_VAR(Int_t, chargedHadronMultiplicity) \
    SIMPLE_VAR(Int_t, chargedMultiplicity) \
    SIMPLE_VAR(Int_t, electronMultiplicity) \
    SIMPLE_VAR(Int_t, muonMultiplicity) \
    SIMPLE_VAR(Int_t, neutralHadronMultiplicity) \
    SIMPLE_VAR(Int_t, neutralMultiplicity) \
    SIMPLE_VAR(Int_t, photonMultiplicity) \
    SIMPLE_VAR(UInt_t, nConstituents) \
    SIMPLE_VAR(Bool_t, passLooseID) \
    SIMPLE_VAR(Bool_t, passTightID) \
    /* Trigger match information */ \
    VECTOR_VAR(std::string, matchedTriggerPaths) \
    /* b-tag discriminators */ \
    B_TAG_DATA() \
    /**/

// see https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/PhysicsTools/PatAlgos/python/tools/jetTools.py#L30
#define B_TAG_DATA() \
    SIMPLE_VAR(Float_t, trackCountingHighEffBJetTags) \
    SIMPLE_VAR(Float_t, trackCountingHighPurBJetTags) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighEffBJetTags) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighPurBJetTags) \
    SIMPLE_VAR(Float_t, jetProbabilityBJetTags) \
    SIMPLE_VAR(Float_t, jetBProbabilityBJetTags) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexBJetTags) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexMVABJetTags) \
    SIMPLE_VAR(Float_t, combinedInclusiveSecondaryVertexBJetTags) \
    SIMPLE_VAR(Float_t, combinedMVABJetTags) \
    /**/


#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Jet, JET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, JetTree, JET_DATA, Jet, "jets", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, JetTree, JET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef JET_DATA

namespace ntuple {
namespace JetID_MVA {
enum Id {
    kTight  = 0,
    kMedium = 1,
    kLoose  = 2
};

inline bool PassId(Int_t puIdBits, Id id) { return ( puIdBits & (1 << id) ) != 0; }
inline bool PassLooseId(Int_t puIdBits) { return PassId(puIdBits, kLoose); }
inline bool PassMediumId(Int_t puIdBits) { return PassId(puIdBits, kMedium); }
inline bool PassTightId(Int_t puIdBits) { return PassId(puIdBits, kTight); }

}
}
