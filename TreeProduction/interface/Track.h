/*! Definiton of ntuple::TrackTree and ntuple::Track classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define TRACK_DATA() \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, etaError) \
    SIMPLE_VAR(Float_t, theta) \
    SIMPLE_VAR(Float_t, thetaError) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, phiError) \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, ptError) \
    SIMPLE_VAR(Float_t, qoverp) \
    SIMPLE_VAR(Float_t, qoverpError) \
    SIMPLE_VAR(Int_t, charge) \
    /* */ \
    SIMPLE_VAR(UInt_t, nValidHits) \
    SIMPLE_VAR(UInt_t, nLostHits) \
    SIMPLE_VAR(Float_t, validFraction) \
    SIMPLE_VAR(Int_t, nValidTrackerHits) \
    SIMPLE_VAR(Int_t, nValidPixelHits) \
    SIMPLE_VAR(Int_t, nValidStripHits) \
    SIMPLE_VAR(Int_t, trackerLayersWithMeasurement) \
    SIMPLE_VAR(Int_t, pixelLayersWithMeasurement) \
    SIMPLE_VAR(Int_t, stripLayersWithMeasurement) \
    /* */ \
    SIMPLE_VAR(Float_t, dxy) \
    SIMPLE_VAR(Float_t, dxyError) \
    SIMPLE_VAR(Float_t, dz) \
    SIMPLE_VAR(Float_t, dzError) \
    SIMPLE_VAR(Float_t, chi2) \
    SIMPLE_VAR(Float_t, ndof) \
    SIMPLE_VAR(Float_t, vx) \
    SIMPLE_VAR(Float_t, vy) \
    SIMPLE_VAR(Float_t, vz) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Track, TRACK_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, TrackTree, TRACK_DATA, Track, "tracks", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, TrackTree, TRACK_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRACK_DATA
