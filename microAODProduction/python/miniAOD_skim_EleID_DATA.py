## @package patTuple
#  Configuration file to produce PAT-tuples and ROOT-tuples for X->HH->bbTauTau analysis.
#
#  \author Claudio Caputo
#
#  Copyright 2015
#
#  This file is part of X->HH->bbTauTau.
#
#  X->HH->bbTauTau is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  X->HH->bbTauTau is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.


import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_reMiniAOD_v1')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


inputSignal = cms.untracked.vstring(
                '/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/6CA1C627-246C-E511-8A6A-02163E014147.root')

process.source = cms.Source("PoolSource",
    fileNames = inputSignal
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

## 'keep recoGsfElectronCores_*_*_PAT',
## 'keep recoSuperClusters_*_*_PAT',
## 'keep patJets_slimmedJetsPuppi__RECO',

skimmedBranches = cms.untracked.vstring(['drop *',
                                         'keep edmTriggerResults_*__HLT',
                                         'keep double_*__RECO',
                                         'keep *_electronMVAValueMapProducer_*_USER', ## Value Map with the electrons ID MVA values
                                         'keep floatedmValueMap_offlineSlimmedPrimaryVertices__RECO',
                                         'keep patPackedTriggerPrescales_patTrigger__RECO',
                                         'keep patElectrons_slimmedElectrons__RECO',
                                         'keep patJets_*__RECO',
                                         'keep patMETs_*__RECO', ## MET
                                         'keep patMuons_slimmedMuons__RECO',
                                         'keep patPackedCandidates_packedPFCandidates__RECO',
                                         'keep patPackedGenParticles_packedGenParticles__RECO',
                                         'keep patTaus_slimmedTaus__RECO',
                                         'keep patTriggerObjectStandAlones_selectedPatTrigger__RECO',
                                         'keep recoGenJets_slimmedGenJets__RECO',
                                         'keep recoGenParticles_prunedGenParticles__RECO',
                                         'keep recoVertexs_offlineSlimmedPrimaryVertices__RECO',
                                         'keep recoVertexCompositePtrCandidates_slimmedSecondaryVertices__RECO',
                                         'keep *_l1extraParticles_*_*',
                                         'keep *_METSignificance_*_USER'])

allBranches = cms.untracked.vstring(['keep *'])

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/cmshome/caputo/HH_bbTauTau/Run2/CMSSW_7_4_12_patch4/src/HHbbTauTau/microAODProduction/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

process.load("RecoMET.METProducers.METSignificance_cfi")
process.load("RecoMET.METProducers.METSignificanceParams_cfi")

process.bbttSkim   = cms.EDFilter("SkimFilterMiniAOD",
                                  vertexSrc  = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  muonSrc  = cms.untracked.InputTag('slimmedMuons'),
                                  electronSrc=cms.untracked.InputTag("slimmedElectrons"),
                                  tauSrc  = cms.untracked.InputTag("slimmedTaus")
                                  )
## Load module for Electron MVA ID
## It will append a value maps the miniAOD, that it's accesible throught a well Handle
## Example code here:
##  https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDwithMVADemo.cc#L99
## process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
##-------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat) ##also compute a maps with the electrons that pass an MVA cut

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
##------------
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("syncTree.root") )

#-------------
# SyncTree Producer
#-------------
process.synctupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc    = cms.InputTag("slimmedTaus"),
                                 muonSrc   = cms.InputTag("slimmedMuons"),
                                 vtxSrc    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 jetSrc    = cms.InputTag("slimmedJets"),
                                 ##pfMETSrc  = cms.InputTag("slimmedMETsNoHF"),
                                 pfMETSrc  = cms.InputTag("slimmedMETs"),
                                 bits      = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects   = cms.InputTag("selectedPatTrigger"),
                                 sampleType = cms.string("Run2015D"),
                                )


process.p = cms.Path(
             #process.electronMVAValueMapProducer*
             process.METSignificance*
             process.egmGsfElectronIDSequence*
             process.bbttSkim*
             process.synctupler
	   	    )
#process.p=cms.Path()

process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('skimBranches_Data.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = skimmedBranches,
    #outputCommands = tausBranch,
    #outputCommands = allBranches,
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.endpath= cms.EndPath(process.out)


