## Configuration file that defines the producer of ROOT-tuple for gen-METs.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

genMETBlock = cms.EDAnalyzer("GenMETBlock",
    verbosity = cms.int32(0),
    genMETSrc = cms.InputTag('genMetTrue')
)
