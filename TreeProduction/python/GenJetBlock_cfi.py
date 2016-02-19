## Configuration file that defines the producer of ROOT-tuple for gen-jets.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

genJetBlock = cms.EDAnalyzer("GenJetBlock",
  verbosity = cms.int32(0),
  genJetSrc = cms.InputTag('ak5GenJets')
)
