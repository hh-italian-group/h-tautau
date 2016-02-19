##  Configuration file that defines a filter to skim events for X->HH->bbTauTau analysis.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

def applySkim(process):
    process.bbttSkim   = cms.EDFilter("SkimFilter",
                                      vertexSrc  = cms.untracked.InputTag('patVertices'),
                                      muonSrc  = cms.untracked.InputTag('slimmedMuons'),
                                      electronSrc=cms.untracked.InputTag("slimmedElectrons"),
                                      tauSrc  = cms.untracked.InputTag("slimmedTaus")
                                      )
    return
