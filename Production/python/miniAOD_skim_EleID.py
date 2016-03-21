# miniAOD skim producer configuration.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


inputSignal = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2CC6C6BD-5303-E511-9F70-0025905964CC.root')
inputTTbar = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root')
inputWJets = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/048FB1EE-33FD-E411-A2BA-0025905A6094.root')
inputDY = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root')
inputSignal_v2 = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/12184969-3DB8-E511-879B-001E67504A65.root')
##inputSignal_v2 = cms.untracked.vstring(
##                '/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/10563B6E-D871-E511-9513-B499BAABD280.root')
## Input files
process.source = cms.Source("PoolSource",
    fileNames = inputSignal_v2
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

## 'keep recoGsfElectronCores_*_*_PAT',
## 'keep recoSuperClusters_*_*_PAT',
## 'keep patJets_slimmedJetsPuppi__PAT',

skimmedBranches = cms.untracked.vstring(['drop *',
                                         'keep edmTriggerResults_*__HLT',
                                         'keep edmTriggerResults_*__PAT',
                                         'keep PileupSummaryInfos_addPileupInfo__HLT',
                                         'keep PileupSummaryInfos_slimmedAddPileupInfo__PAT',
                                         'keep double_*__RECO',
                                         'keep recoBeamSpot_offlineBeamSpot__RECO',
                                         'keep *_electronMVAValueMapProducer_*_USER', ## Value Map with the electrons ID MVA values
                                         'keep floatedmValueMap_offlineSlimmedPrimaryVertices__PAT',
                                         'keep patPackedTriggerPrescales_patTrigger__PAT',
                                         'keep patElectrons_slimmedElectrons__PAT',
                                         'keep patJets_*__PAT',
                                         'keep patMETs_*__PAT', ## MET
                                         'keep patMuons_slimmedMuons__PAT',
                                         'keep patPackedCandidates_packedPFCandidates__PAT',
                                         'keep patPackedGenParticles_packedGenParticles__PAT',
                                         'keep patTaus_slimmedTaus__PAT',
                                         'keep patTriggerObjectStandAlones_selectedPatTrigger__PAT',
                                         'keep recoGenJets_slimmedGenJets__PAT',
                                         'keep recoGenParticles_prunedGenParticles__PAT',
                                         'keep recoVertexs_offlineSlimmedPrimaryVertices__PAT',
                                         'keep recoVertexCompositePtrCandidates_slimmedSecondaryVertices__PAT',
                                         'keep *_l1extraParticles_*_*',
                                         'keep *_METSignificance_*_USER'])

allBranches = cms.untracked.vstring(['keep *'])

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

process.p = cms.Path(
             #process.electronMVAValueMapProducer*
             process.bbttSkim*
             process.METSignificance*
             process.egmGsfElectronIDSequence

	   	    )
#process.p=cms.Path()

process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('skimBranches_Signal_miniV2.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = skimmedBranches,
    #outputCommands = tausBranch,
    #outputCommands = allBranches,
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.endpath= cms.EndPath(process.out)


