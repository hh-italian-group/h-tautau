##  Configuration file that defines the producer of ROOT-tuple for vertices.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.int32(1),
  vertexSrc = cms.InputTag('patVertices')
)
