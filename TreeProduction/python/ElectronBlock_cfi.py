## Configuration file that defines the producer of ROOT-tuple for electrons.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity       = cms.int32(0),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  trackSrc        = cms.InputTag('generalTracks'),
  vertexSrc       = cms.InputTag('patVertices'),
  electronSrc     = cms.InputTag('patElectronsWithEmbeddedVariables')
)
