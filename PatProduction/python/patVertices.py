##  Configuration file that defines parameters related to PAT Vertex objects.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

def applyVertexParameters(process):
    process.patVertices = cms.EDProducer('PatVertexProducer',
		inputTag = cms.InputTag('offlinePrimaryVertices')
    )

    return
