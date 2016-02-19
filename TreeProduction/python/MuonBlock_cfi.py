## Configuration file that defines the producer of ROOT-tuple for muons.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  muonSrc         = cms.InputTag('patMuonsWithEmbeddedVariables'),
  vertexSrc       = cms.InputTag('patVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  beamSpotCorr    = cms.bool(True),
  muonID          = cms.string('GlobalMuonPromptTight')
)
