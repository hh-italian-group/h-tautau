/*! Definiton of ntuple::ElectronTree and ntuple::Electron classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define ELECTRON_DATA() \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, mass) \
    SIMPLE_VAR(Float_t, caloEnergy) \
    SIMPLE_VAR(Float_t, caloEnergyError) \
    /* Origin */ \
    SIMPLE_VAR(Bool_t, ecalDriven) \
    SIMPLE_VAR(Bool_t, hasGsfTrack) \
    /* Gsf based electrons */ \
    SIMPLE_VAR(Float_t, trackPt) \
    SIMPLE_VAR(Float_t, trackPtError) \
    SIMPLE_VAR(Int_t, charge) \
    SIMPLE_VAR(Int_t, pixHits) \
    SIMPLE_VAR(Int_t, trkHits) \
    SIMPLE_VAR(UInt_t, nValidHits) \
    SIMPLE_VAR(Float_t, trkD0) \
    SIMPLE_VAR(Float_t, trkD0Error) \
    SIMPLE_VAR(Int_t, missingHits) \
    /* ID variables */ \
    SIMPLE_VAR(Float_t, hcalOverEcal) \
    SIMPLE_VAR(Float_t, hcalDepth1OverEcal) \
    SIMPLE_VAR(Float_t, eSuperClusterOverP) \
    SIMPLE_VAR(Float_t, sigmaEtaEta) \
    SIMPLE_VAR(Float_t, sigmaIEtaIEta) \
    SIMPLE_VAR(Float_t, deltaPhiTrkSC) \
    SIMPLE_VAR(Float_t, deltaEtaTrkSC) \
    SIMPLE_VAR(Int_t, classification) \
    SIMPLE_VAR(Float_t, e1x5overe5x5) \
    SIMPLE_VAR(Float_t, e2x5overe5x5) \
    /* Iso variables */ \
    SIMPLE_VAR(Float_t, isoEcal03) \
    SIMPLE_VAR(Float_t, isoHcal03) \
    SIMPLE_VAR(Float_t, isoTrk03) \
    SIMPLE_VAR(Float_t, isoEcal04) \
    SIMPLE_VAR(Float_t, isoHcal04) \
    SIMPLE_VAR(Float_t, isoTrk04) \
    SIMPLE_VAR(Float_t, isoRel03) \
    SIMPLE_VAR(Float_t, isoRel04) \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vx) \
    SIMPLE_VAR(Float_t, vy) \
    SIMPLE_VAR(Float_t, vz) \
    /* SC associated with electron */ \
    SIMPLE_VAR(Float_t, scEn) \
    SIMPLE_VAR(Float_t, scEta) \
    SIMPLE_VAR(Float_t, scPhi) \
    SIMPLE_VAR(Float_t, scET) \
    SIMPLE_VAR(Float_t, scRawEnergy) \
    /* Vertex association variables */ \
    SIMPLE_VAR(Float_t, vtxDist3D) \
    SIMPLE_VAR(Int_t, vtxIndex) \
    SIMPLE_VAR(Float_t, vtxDistZ) \
    SIMPLE_VAR(Float_t, relIso) \
    SIMPLE_VAR(Float_t, pfRelIso) \
    /* PFlow isolation variable */ \
    SIMPLE_VAR(Float_t, chargedHadronIso) \
    SIMPLE_VAR(Float_t, neutralHadronIso) \
    SIMPLE_VAR(Float_t, photonIso) \
    /* IP information */ \
    SIMPLE_VAR(Float_t, dB) \
    SIMPLE_VAR(Float_t, edB) \
    SIMPLE_VAR(Float_t, dB3d) \
    SIMPLE_VAR(Float_t, edB3d) \
    SIMPLE_VAR(Int_t, nBrems) \
    SIMPLE_VAR(Float_t, fbrem) \
    SIMPLE_VAR(Float_t, dist_vec) \
    SIMPLE_VAR(Float_t, dCotTheta) \
    SIMPLE_VAR(Bool_t, hasMatchedConversion) \
    /* MVA */ \
    SIMPLE_VAR(Float_t, mva) \
    SIMPLE_VAR(Float_t, mvaPOGTrig) \
    SIMPLE_VAR(Float_t, mvaPOGNonTrig) \
    SIMPLE_VAR(Bool_t, mvaPreselection) \
    SIMPLE_VAR(Bool_t, isTriggerElectron) \
    SIMPLE_VAR(Float_t, isoMVA) \
    SIMPLE_VAR(Float_t, pfRelIso03v1) \
    SIMPLE_VAR(Float_t, pfRelIso03v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v1) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v3) \
    SIMPLE_VAR(Float_t, pfRelIso04v1) \
    SIMPLE_VAR(Float_t, pfRelIso04v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v1) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v3) \
    SIMPLE_VAR(Float_t, pfRelIso03) \
    SIMPLE_VAR(Float_t, pfRelIso04) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04) \
    SIMPLE_VAR(UInt_t, fidFlag) \
    /* Trigger match information */ \
    VECTOR_VAR(std::string, matchedTriggerPaths)
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Electron, ELECTRON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, ElectronTree, ELECTRON_DATA, Electron, "electrons", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, ElectronTree, ELECTRON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef ELECTRON_DATA
