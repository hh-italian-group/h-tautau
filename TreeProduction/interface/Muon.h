/*! Definiton of ntuple::MuonTree and ntuple::Muon classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define MUON_DATA() \
    /* Muon ID */ \
    SIMPLE_VAR(Bool_t, isTrackerMuon) \
    SIMPLE_VAR(Bool_t, isPFMuon) \
    SIMPLE_VAR(Bool_t, isGlobalMuonPromptTight) \
    SIMPLE_VAR(Bool_t, isTightMuon) \
    SIMPLE_VAR(Bool_t, isHighPtMuon) \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, mass) \
    SIMPLE_VAR(Float_t, ptError) \
    SIMPLE_VAR(Int_t, charge) \
    /* Track info */ \
    SIMPLE_VAR(Float_t, trkD0) \
    SIMPLE_VAR(Float_t, trkD0Error) \
    SIMPLE_VAR(Float_t, trkDz) \
    SIMPLE_VAR(Float_t, trkDzError) \
    SIMPLE_VAR(Float_t, globalChi2) \
    SIMPLE_VAR(Bool_t, passID) \
    SIMPLE_VAR(Int_t, pixHits) \
    SIMPLE_VAR(Int_t, trkHits) \
    SIMPLE_VAR(Int_t, muonHits) \
    SIMPLE_VAR(Int_t, matches) \
    SIMPLE_VAR(Int_t, trackerLayersWithMeasurement) \
    /* Isolation */ \
    SIMPLE_VAR(Float_t, trkIso) \
    SIMPLE_VAR(Float_t, ecalIso) \
    SIMPLE_VAR(Float_t, hcalIso) \
    SIMPLE_VAR(Float_t, hoIso) \
    SIMPLE_VAR(Float_t, relIso) \
    SIMPLE_VAR(Float_t, pfRelIso) \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vtxDist3D) \
    SIMPLE_VAR(Int_t, vtxIndex) \
    SIMPLE_VAR(Float_t, vtxDistZ) \
    SIMPLE_VAR(Float_t, vx) \
    SIMPLE_VAR(Float_t, vy) \
    SIMPLE_VAR(Float_t, vz)\
    /* PV2D*/ \
    SIMPLE_VAR(Float_t, dB) \
    SIMPLE_VAR(Float_t, edB) \
    /* PV3D */ \
    SIMPLE_VAR(Float_t, dB3d) \
    SIMPLE_VAR(Float_t, edB3d) \
    /* UW Recommendation*/ \
    SIMPLE_VAR(Bool_t, isAllArbitrated) \
    SIMPLE_VAR(Int_t, nChambers) \
    SIMPLE_VAR(Int_t, nMatches) \
    SIMPLE_VAR(Int_t, nMatchedStations) \
    SIMPLE_VAR(UInt_t, stationMask) \
    SIMPLE_VAR(UInt_t, stationGapMaskDistance)\
    SIMPLE_VAR(UInt_t, stationGapMaskPull) \
    SIMPLE_VAR(Int_t, muonID) \
    /* MVA */ \
    SIMPLE_VAR(Float_t, idMVA) \
    SIMPLE_VAR(Float_t, isoRingsMVA) \
    SIMPLE_VAR(Float_t, isoRingsRadMVA) \
    SIMPLE_VAR(Float_t, idIsoCombMVA) \
    /* Iso variables*/ \
    SIMPLE_VAR(Float_t, pfRelIso03v1) \
    SIMPLE_VAR(Float_t, pfRelIso03v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v1) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v2) \
    SIMPLE_VAR(Float_t, pfRelIso04v1) \
    SIMPLE_VAR(Float_t, pfRelIso04v2) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v1) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v2) \
    /* Trigger match information */ \
    VECTOR_VAR(std::string, matchedTriggerPaths)
    /**/


#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Muon, MUON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, MuonTree, MUON_DATA, Muon, "muons", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, MuonTree, MUON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef MUON_DATA
