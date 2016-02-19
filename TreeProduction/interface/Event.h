/*! Definiton of ntuple::EventTree and ntuple::Event classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define EVENT_DATA() \
    SIMPLE_VAR(UInt_t, run) \
    SIMPLE_VAR(UInt_t, EventId) \
    SIMPLE_VAR(UInt_t, lumis) \
    SIMPLE_VAR(Int_t, bunch) \
    SIMPLE_VAR(Int_t, orbit) \
    SIMPLE_VAR(UInt_t, unixTime) \
    SIMPLE_VAR(UInt_t, microsecondOffset) \
    SIMPLE_VAR(Bool_t, isdata) \
    SIMPLE_VAR(Bool_t, isPhysDeclared)  \
    SIMPLE_VAR(Bool_t, isBPTX0) \
    SIMPLE_VAR(Bool_t, isBSCMinBias) \
    SIMPLE_VAR(Bool_t, isBSCBeamHalo) \
    SIMPLE_VAR(Bool_t, isPrimaryVertex) \
    VECTOR_VAR(Int_t, nPU) \
    VECTOR_VAR(Int_t, bunchCrossing) \
    VECTOR_VAR(Float_t, trueNInt) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Event, EVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, EventTree, EVENT_DATA, Event, "events", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, EventTree, EVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef EVENT_DATA
