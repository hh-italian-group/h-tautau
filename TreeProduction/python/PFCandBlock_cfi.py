## Configuration file that defines the producer of ROOT-tuple for PF candidates.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

pfCandBlock = cms.EDAnalyzer("PFCandBlock",
    srcPFCandidates = cms.InputTag('particleFlow'),
)
