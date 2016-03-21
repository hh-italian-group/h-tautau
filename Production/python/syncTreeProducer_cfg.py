# Sync tree producer configuration.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauSyncTree2015")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globalTag',
                  'MCRUN2_74_V9::All',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Global Tag to use. Default: MCRUN2_74_V9::All")
options.register ('fileList',
                  'fileList.txt',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")
options.register ('fileNamePrefix',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Prefix to add to input file names.")
options.register ('includeSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")
options.register ('isData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")
options.register ('runOnCrab', 
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Indicates if script will be executed on CRAB.")
options.register ('sampleType', 
                  'Fall15MC',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
options.register ('computeHT',
                'False',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                "Compute HT variable and HT binning")


options.parseArguments()

#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------

from HHbbTauTau.RunTools.readFileList import *

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring())

if not options.runOnCrab:
  readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
  process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.Geometry.GeometryIdeal_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag and Lumi
#-------------
import FWCore.PythonUtilities.LumiList as LumiList
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'

runOnData = options.isData

if runOnData:
  process.GlobalTag.globaltag = '76X_dataRun2_v15'
  process.source.lumisToProcess = LumiList.LumiList(filename = '/cmshome/caputo/HH_bbTauTau/Run2/CMSSW_7_4_12_patch4/src/HHbbTauTau/microAODProduction/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()
else:
  process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile) )

#-------------
# SyncTree Producer
#-------------
#from HHbbTauTau.microAODProduction.syncNtupler_cfi import syncNtupler

process.syncNtupler = cms.EDAnalyzer('SyncTreeProducer',
                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc     = cms.InputTag("slimmedTaus"),
                                 muonSrc    = cms.InputTag("slimmedMuons"),
                                 vtxSrc     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 jetSrc     = cms.InputTag("slimmedJets"),
                                 PUInfo    = cms.InputTag("slimmedAddPileupInfo"),
                                 ##pfMETSrc  = cms.InputTag("slimmedMETsNoHF"),
                                 pfMETSrc   = cms.InputTag("slimmedMETs"),
                                 bits       = cms.InputTag("TriggerResults","","HLT"),
                                 prescales  = cms.InputTag("patTrigger"),
                                 objects    = cms.InputTag("selectedPatTrigger"),
                                 metCov     = cms.InputTag("METSignificance","METCovariance"),
                                 sampleType = cms.string(options.sampleType),
				)

process.p = cms.Path(process.synctupler)
