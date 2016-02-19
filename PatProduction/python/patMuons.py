##  Configuration file that defines parameters related to PAT Muon objects.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso
from CommonTools.ParticleFlow.pfParticleSelection_cff import pfParticleSelectionSequence

def applyMuonParameters(process):
    process.muIsoSequence       = setupPFMuonIso(process,'muons')
    process.pfParticleSelectionSequence = pfParticleSelectionSequence

    process.patMuons.isoDeposits = cms.PSet(
        pfAllParticles   = cms.InputTag("muPFIsoDepositPUPFIso"),              # all PU   CH+MU+E
        pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"),         # all noPU CH
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"),         # all NH
        pfPhotons        = cms.InputTag("muPFIsoDepositGammaPFIso"),           # all PH
        user = cms.VInputTag( cms.InputTag("muPFIsoDepositChargedAllPFIso") )  # all noPU CH+MU+E
        )

    process.patMuons.isolationValues = cms.PSet(
        pfAllParticles   = cms.InputTag("muPFIsoValuePU04PFIso"),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
        pfPhotons        = cms.InputTag("muPFIsoValueGamma04PFIso"),
        user = cms.VInputTag( cms.InputTag("muPFIsoValueChargedAll04PFIso") )
        )

    process.patMuons.embedTrack = cms.bool(True)

    process.patMuonsWithEmbeddedVariables = cms.EDProducer('MuonsUserEmbedded',
        muonTag = cms.InputTag("patMuonsTriggerMatch"),
        vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS")
        )

    return
