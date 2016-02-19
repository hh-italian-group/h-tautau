## Configuration file that defines the producer of ROOT-tuple for gen-particles.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

genEventBlock = cms.EDAnalyzer("GenEventBlock",
    genEventSrc = cms.InputTag('generator','minVisPtFilter','EmbeddedRECO'),
    lheProductSrc = cms.InputTag('source')
)
