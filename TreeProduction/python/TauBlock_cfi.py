## Configuration file that defines the producer of ROOT-tuple for taus.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
    patTauSrc = cms.InputTag('patTausTriggerMatch'),
)
