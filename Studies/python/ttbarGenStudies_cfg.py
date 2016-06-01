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
options = VarParsing('analysis')
options.register ('isData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")
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

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

runOnData = options.isData

# Latest JEC
if runOnData:
  process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
  isMC = False
  #process.source.lumisToProcess = LumiList.LumiList(filename = '../json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()
else:
  process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
  isMC = True


## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

inputSignal_v2 = cms.untracked.vstring("file:768F5AFB-D771-E511-9ABD-B499BAABD280.root")

DYSample = cms.untracked.vstring("/store/user/ccaputo/HHbbtautau/Run2/DYSample_forHT.root")

TTBar = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/0007DBD0-2ED2-E511-AD0D-20CF3019DEF5.root')
SyncSignal = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/B2FF8F77-3DB8-E511-B743-001E6757F1D4.root')
Radion300  = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/821F9C9F-28B8-E511-93CF-003048D2BF3E.root')
SyncSignal_fewEvents = cms.untracked.vstring(("root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/E0B9088F-3DB8-E511-AFFD-001EC9ADCD52.root",
											  "root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/EAD78CAB-33B8-E511-8E8E-20CF3027A5BF.root"))
##inputSignal_v2 = cms.untracked.vstring(
##                '/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/10563B6E-D871-E511-9513-B499BAABD280.root')
## Input files
process.source = cms.Source("PoolSource",
    fileNames = Radion300#,
                                 #eventsToProcess = cms.untracked.VEventRange('1:449465','1:475952')
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# from h-tautau.Production.skimmedBranches_cff import *
#
# if options.computeHT and not options.isData:
#     skimmedBranches = cms.untracked.vstring(BaseMCBranches+
#                                             ['keep LHEEventProduct_externalLHEProducer__LHE'])
# if not options.computeHT and not options.isData:
#     skimmedBranches = cms.untracked.vstring(BaseMCBranches)
#
# if options.isData:
#     skimmedBranches = cms.untracked.vstring(BaseDATABranches)

ttbarStudiesBrances = cms.untracked.vstring(['drop *',
                                            'keep patPackedCandidates_packedPFCandidates__PAT',
                                            'keep patPackedGenParticles_packedGenParticles__PAT',
                                            'keep recoGenJets_slimmedGenJets__PAT',
                                            'keep recoGenParticles_prunedGenParticles__PAT',
                                             'keep GenEventInfoProduct_generator__SIM',
                                             'keep patJets_*__PAT',
                                            'keep patMETs_*__PAT']
                                             )

DropAllBranches = cms.untracked.vstring(['drop *'])

process.load("RecoMET.METProducers.METSignificance_cfi")
process.load("RecoMET.METProducers.METSignificanceParams_cfi")

process.bbttSkim   = cms.EDFilter("SkimFilterMiniAOD",
                                  vertexSrc  = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  muonSrc  = cms.untracked.InputTag('slimmedMuons'),
                                  electronSrc=cms.untracked.InputTag("slimmedElectrons"),
                                  tauSrc  = cms.untracked.InputTag("slimmedTaus")
                                  )
process.ttbarAnalyzer   = cms.EDAnalyzer("TTBarGenAnalyzer",
                                  pruned  = cms.InputTag('prunedGenParticles'),
                                  packed  = cms.InputTag('packedGenParticles'),
                                  genJet  = cms.InputTag('slimmedGenJets'),
                                  isSignal = cms.bool(True),
                                  )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                    maxEventsToPrint = cms.untracked.int32(1),
                                    printVertex = cms.untracked.bool(False),
                                    printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
                                    src = cms.InputTag("prunedGenParticles")
                                    )
##------------
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root") )

process.p = cms.Path(
             # process.METSignificance*
#              process.egmGsfElectronIDSequence*
#              process.electronMVAValueMapProducer*
             process.bbttSkim*
             process.printTree*
             process.ttbarAnalyzer
             #(process.syncNtupler_mutau + process.syncNtupler_etau)
	   	    )


process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('microAOD.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #outputCommands = skimmedBranches,
    #outputCommands = HTBinBranches,
    outputCommands = ttbarStudiesBrances,
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.endpath= cms.EndPath(process.out)
