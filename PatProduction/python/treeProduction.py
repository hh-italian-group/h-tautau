##  Configuration file that defines the sequence to produce ROOT-tuples for X->HH->bbTauTau analysis.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

def addTreeSequence(process, includeSim, treeOutput):
    process.TFileService = cms.Service("TFileService", fileName = cms.string(treeOutput) )
    process.load("HHbbTauTau.TreeProduction.TreeContentConfig_cff")

    process.mainTreeContentSequence = cms.Sequence(
        process.eventBlock
      + process.vertexBlock
      + process.electronBlock
      + process.jetBlock
      + process.metBlock
      + process.muonBlock
      + process.tauBlock
      + process.pfCandBlock
      + process.triggerBlock
      + process.triggerObjectBlock
    )

    process.simTreeContentSequence = cms.Sequence()

    if includeSim:
        process.simTreeContentSequence = cms.Sequence(process.genParticleBlock + process.genEventBlock + process.genMETBlock)

    return
