## Configuration file that defines the producer of ROOT-tuple for events.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

eventBlock = cms.EDAnalyzer("EventBlock",
    verbosity = cms.int32(0),
    l1InputTag  = cms.InputTag('gtDigis'),
    vertexInputTag = cms.InputTag('patVertices'),
    vertexMinimumNDOF = cms.uint32(4),
    vertexMaxAbsZ = cms.double(24.),
    vertexMaxd0 = cms.double(2.),
    hpTrackThreshold = cms.double(0.25)
)
