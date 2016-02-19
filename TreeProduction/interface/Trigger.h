/*! Definiton of ntuple::TriggerTree and ntuple::Trigger classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define TRIGGER_DATA() \
    VECTOR_VAR(Bool_t, l1physbits) \
    VECTOR_VAR(Bool_t, l1techbits) \
    VECTOR_VAR(std::string, hltpaths) \
    VECTOR_VAR(Bool_t, hltresults) \
    VECTOR_VAR(UInt_t, hltprescales) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Trigger, TRIGGER_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, TriggerTree, TRIGGER_DATA, Trigger, "triggers", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, TriggerTree, TRIGGER_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRIGGER_DATA
