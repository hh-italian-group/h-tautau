/*! Definiton of ntuple::GenJetTree and ntuple::GenJet classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define GENJET_DATA() \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, mass) \
    SIMPLE_VAR(Float_t, emEnergy) \
    SIMPLE_VAR(Float_t, hadEnergy) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, GenJet, GENJET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, GenJetTree, GENJET_DATA, GenJet, "genJets", true)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, GenJetTree, GENJET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GENJET_DATA
