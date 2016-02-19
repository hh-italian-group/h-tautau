/*! Definiton of ntuple::GenParticleTree and ntuple::GenParticle classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define GEN_PARTICLE_DATA() \
    SIMPLE_VAR(Int_t, PdgId) \
    SIMPLE_VAR(Int_t, Status) \
    SIMPLE_VAR(Int_t, Charge) \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, mass) \
    SIMPLE_VAR(Float_t, X) \
    SIMPLE_VAR(Float_t, Y) \
    SIMPLE_VAR(Float_t, Z) \
    VECTOR_VAR(UInt_t, Mother_Indexes) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, GenParticle, GEN_PARTICLE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, GenParticleTree, GEN_PARTICLE_DATA, GenParticle, "genParticles", true)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, GenParticleTree, GEN_PARTICLE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GEN_PARTICLE_DATA
