## Configuration file that defines the producer of ROOT-tuple for jets.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

jetBlock = cms.EDAnalyzer("JetBlock",
    jetSrc = cms.InputTag('PatJetsWithEmbeddedVariables')
)
