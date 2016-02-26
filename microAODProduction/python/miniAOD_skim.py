## @package patTuple
#  Configuration file to produce PAT-tuples and ROOT-tuples for X->HH->bbTauTau analysis.
#
#  \author
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
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )


inputSignal = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2CC6C6BD-5303-E511-9F70-0025905964CC.root')
inputTTbar = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root')
inputWJets = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/048FB1EE-33FD-E411-A2BA-0025905A6094.root')
inputDY = cms.untracked.vstring(
                '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root')

AZh_minAODv2 = cms.untracked.vstring(
                '/store/mc/RunIISpring15MiniAODv2/AToZhToLLTauTau_M-220_13TeV_madgraph_4f_LO/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/E47E8231-746F-E511-93EF-34E6D7BDDEDB.root')
inputSignal_v2 = cms.untracked.vstring(
                '/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/18442D74-D671-E511-99AF-00266CFAE748.root')

## Input files
process.source = cms.Source("PoolSource",
    fileNames = AZh_minAODv2
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
                                         'keep double_*__RECO',
                                         'keep recoBeamSpot_offlineBeamSpot__RECO',
                                         'keep floatedmValueMap_offlineSlimmedPrimaryVertices__PAT',
                                         'keep patPackedTriggerPrescales_patTrigger__PAT',
                                         'keep patElectrons_slimmedElectrons__PAT',
                                         'keep patJets_*__PAT',
                                         'keep patMETs_*__PAT',
                                         'keep patMuons_slimmedMuons__PAT',
                                         'keep patPackedCandidates_packedPFCandidates__PAT',
                                         'keep patPackedGenParticles_packedGenParticles__PAT',
                                         'keep patTaus_slimmedTaus__PAT',
                                         'keep patTriggerObjectStandAlones_selectedPatTrigger__PAT',
                                         'keep recoGenJets_slimmedGenJets__PAT',
                                         'keep recoGenParticles_prunedGenParticles__PAT',
                                         'keep recoVertexs_offlineSlimmedPrimaryVertices__PAT',
                                         'keep recoVertexCompositePtrCandidates_slimmedSecondaryVertices__PAT',
                                         'keep *_l1extraParticles_*_*'])

allBranches = cms.untracked.vstring(['keep *'])

tausBranch = cms.untracked.vstring(['keep patTaus_slimmedTaus__PAT'])



process.bbttSkim   = cms.EDFilter("SkimFilterMiniAOD",
                                  vertexSrc  = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  muonSrc  = cms.untracked.InputTag('slimmedMuons'),
                                  electronSrc=cms.untracked.InputTag("slimmedElectrons"),
                                  tauSrc  = cms.untracked.InputTag("slimmedTaus")
                                  )

process.p = cms.Path(process.bbttSkim)
#process.p=cms.Path()

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('allBranches_miniAOD_AZh.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #outputCommands = skimmedBranches
    #outputCommands = tausBranch
    outputCommands = allBranches
)
process.endpath= cms.EndPath(process.out)


