# Produce EventTuple for all channels.
# Based on https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/HepMCCandAlgos/test/printTree.py
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('inputFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Input root files to process.")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of events to process.")
options.parseArguments()

process = cms.Process("PrintGenTruth")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source( "PoolSource", fileNames = cms.untracked.vstring(options.inputFile))
if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

process.printTree = cms.EDAnalyzer("PrintGenTruth",
    genParticles = cms.InputTag("prunedGenParticles"),
    lheEventProduct = cms.InputTag('externalLHEProducer')
)

process.p = cms.Path(process.printTree)
