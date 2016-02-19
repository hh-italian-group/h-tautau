##  Configuration file that defines command-line options for PATtoople production for X->HH->bbTauTau analysis.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from HHbbTauTau.RunTools.readFileList import *

import sys

options = VarParsing('analysis')

def parseAndApplyOptions(process) :
    options.register ('globalTag', 'START53_V7A::All', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Global Tag to use.")
    options.register ('isMC', True, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Sample Type: MC or data.")
    options.register ('isEmbedded', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Indicates if sample is made using embedding technique.")
    options.register ('runOnCrab', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Indicates if script will be executed on CRAB.")
    options.register ('fileList', 'fileList.txt', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "List of root files to process.")
    options.register ('fileNamePrefix', '', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Prefix to add to input file names.")
    options.register ('includeSim', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Include Sim information.")
    options.register ('keepPat', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Keep PAT information in the output file.")
    options.register ('runTree', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Run TREE production sequence.")
    options.register ('treeOutput', 'Tree.root', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Tree root file.")

    if len(sys.argv) > 0:
        last = sys.argv.pop()
        sys.argv.extend(last.split(","))

    options.parseArguments()

    process.GlobalTag.globaltag = options.globalTag

    if not options.runOnCrab:
        readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
        ## Output file
#        from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
#        process.OUT = cms.OutputModule("PoolOutputModule",
#            fileName = cms.untracked.string('test.root'),
#            outputCommands = cms.untracked.vstring(['keep *'])
#        )
#        process.endpath= cms.EndPath(process.OUT)
        process.maxEvents.input = options.maxEvents

    return
