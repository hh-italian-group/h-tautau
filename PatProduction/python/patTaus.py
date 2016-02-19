##  Configuration file that defines parameters related to PAT Tau objects.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.tauTools import *

## Should work for both standard and high pt taus
def applyTauParameters(process):
    # Include tau reconstruction sequence.
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    switchToPFTauHPS(process)

    process.patTaus.embedLeadPFCand = cms.bool(True)
    process.patTaus.embedSignalPFCands = cms.bool(True)
    process.patTaus.embedIsolationPFCands = cms.bool(True)
    process.patTaus.embedLeadTrack = cms.bool(True)
    process.patTaus.embedSignalTracks = cms.bool(True)
    process.patTaus.embedIsolationTracks = cms.bool(True)
    process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFGammaCands = cms.bool(True)
    process.patTaus.embedSignalPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFGammaCands = cms.bool(True)
    process.patTaus.embedLeadPFChargedHadrCand = cms.bool(True)
    process.patTaus.embedLeadPFNeutralCand = cms.bool(True)

    return
