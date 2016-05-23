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
options.register ('ReRunJEC',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Re-Run Jet Energy Corrections. Default: False")                  
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

process = cms.Process("TreeProducer")

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

Diboson = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/1A826380-CB6D-E511-BCFA-0025901D4D6E.root')
SyncSignal = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/12184969-3DB8-E511-879B-001E67504A65.root')

SyncSignal_fewEvents = cms.untracked.vstring(("root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/E0B9088F-3DB8-E511-AFFD-001EC9ADCD52.root",
											  "root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/EAD78CAB-33B8-E511-8E8E-20CF3027A5BF.root"))
##inputSignal_v2 = cms.untracked.vstring(
##                '/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/10563B6E-D871-E511-9513-B499BAABD280.root')
## Input files
process.source = cms.Source("PoolSource",
    fileNames = SyncSignal#,
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

allBranches = cms.untracked.vstring(['drop *'])

### ReRun JEC ===========================================================================================
if options.ReRunJEC:
   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
   process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
   src = cms.InputTag("slimmedJets"),
   levels = ['L1FastJet', 
           'L2Relative', 
           'L3Absolute'],
   payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

   if not isMC:
      process.patJetCorrFactorsReapplyJEC.levels.append("L2L3Residual")

   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
   process.patJetsReapplyJEC = patJetsUpdated.clone(
   jetSource = cms.InputTag("slimmedJets"),
   jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
   )
   process.JECsequence = cms.Sequence (process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC)
   JetCollection = cms.InputTag("patJetsReapplyJEC::TreeProducer")
else:
   process.JECsequence = cms.Sequence ()
   JetCollection = cms.InputTag("slimmedJets")      

### END ReRun JEC ======================================================================================


process.load("RecoMET.METProducers.METSignificance_cfi")
process.load("RecoMET.METProducers.METSignificanceParams_cfi")

if options.ReRunJEC:
   process.METCorrSignificance = process.METSignificance.clone(
      srcPfJets = cms.InputTag('patJetsReapplyJEC::TreeProducer'),
      srcMet = cms.InputTag('patpfMETT1::TreeProducer')
      )



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

process.syncNtupler_etau = cms.EDAnalyzer('SyncTreeProducer_etau',
                                            genParticles  = cms.InputTag("genParticles"),
                                            electronSrc   = cms.InputTag("slimmedElectrons"),
                                            eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                            eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                            eleCutBasedVeto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
					    tauSrc        = cms.InputTag("slimmedTaus"),
                                            muonSrc       = cms.InputTag("slimmedMuons"),
                                            vtxSrc        = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            jetSrc        = JetCollection,
                                            PUInfo    = cms.InputTag("slimmedAddPileupInfo"),
                                            ##pfMETSrc       = cms.InputTag("slimmedMETsNoHF"),
                                            pfMETSrc         = cms.InputTag("slimmedMETs"),
                                            bits             = cms.InputTag("TriggerResults","","HLT"),
                                            prescales        = cms.InputTag("patTrigger"),
                                            objects          = cms.InputTag("selectedPatTrigger"),
                                            metCov     = cms.InputTag("METSignificance","METCovariance"),
                                            lheEventProducts = cms.InputTag("externalLHEProducer"),
                                            genEventInfoProduct = cms.InputTag("generator"),
                                            pruned  = cms.InputTag('prunedGenParticles'),
                                            l1JetParticleProduct = cms.InputTag("l1extraParticles","IsoTau"),
                                            isMC             = cms.bool(isMC),
                                            HTBinning        = cms.bool(options.computeHT),
                                            sampleType = cms.string(options.sampleType),
                                            )
                                           

process.p = cms.Path(
             process.JECsequence *
             process.METSignificance*
             process.egmGsfElectronIDSequence*
             process.electronMVAValueMapProducer*
             process.bbttSkim*
             process.syncNtupler_etau
	   	    )


process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('microAOD.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #outputCommands = skimmedBranches,
    #outputCommands = HTBinBranches,
    outputCommands = allBranches,
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.endpath= cms.EndPath(process.out)
