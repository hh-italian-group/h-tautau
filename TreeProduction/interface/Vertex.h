/*! Definiton of ntuple::VertexTree and ntuple::Vertex classes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

#define VERTEX_DATA() \
    SIMPLE_VAR(Float_t, x) \
    SIMPLE_VAR(Float_t, y) \
    SIMPLE_VAR(Float_t, z) \
    SIMPLE_VAR(Float_t, xErr) \
    SIMPLE_VAR(Float_t, yErr) \
    SIMPLE_VAR(Float_t, zErr) \
    SIMPLE_VAR(Float_t, rho) \
    SIMPLE_VAR(Float_t, chi2) \
    SIMPLE_VAR(Float_t, ndf) \
    SIMPLE_VAR(UInt_t, ntracks) \
    SIMPLE_VAR(UInt_t, ntracksw05) \
    SIMPLE_VAR(Bool_t, isfake) \
    SIMPLE_VAR(Bool_t, isvalid) \
    SIMPLE_VAR(Float_t, sumPt) \
    SIMPLE_VAR(Float_t, sumPtSquared) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Vertex, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, VertexTree, VERTEX_DATA, Vertex, "vertices", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, VertexTree, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef VERTEX_DATA
