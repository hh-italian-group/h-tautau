##  Configuration file that defines the producer of ROOT-tuple for general tracks.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

trackBlock = cms.EDAnalyzer("TrackBlock",
  verbosity       = cms.int32(0),
  trackSrc        = cms.InputTag('generalTracks'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot')
)
