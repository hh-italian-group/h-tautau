# Sync tree producer.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

process.syncNtupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc           = cms.InputTag("slimmedTaus"),
                                 muonSrc          = cms.InputTag("slimmedMuons"),
                                 vtxSrc           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 jetSrc           = cms.InputTag("slimmedJets"),
                                 ##pfMETSrc       = cms.InputTag("slimmedMETsNoHF"),
                                 pfMETSrc         = cms.InputTag("slimmedMETs"),
                                 bits             = cms.InputTag("TriggerResults","","HLT"),
                                 prescales        = cms.InputTag("patTrigger"),
                                 objects          = cms.InputTag("selectedPatTrigger"),
                                 lheEventProducts = cms.InputTag("externalLHEProducer"),
                                 HTBinning        = cms.bool(options.computeHT),
                                 sampleType = cms.string(options.sampleType),
                                )
