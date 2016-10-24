# Produce EventTuple for all channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('globalTag', '76X_mcRun2_asymptotic_RunIIFall15DR76_v1', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Global Tag to use.")
options.register('sampleType', 'Fall15MC', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
options.register('ReRunJEC', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Re-run Jet Energy Corrections. Default: False")
options.register('applyTriggerMatch', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Apply trigger matching for signal objects. Default: True")
options.register('fileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                  "Prefix to add to input file names.")
options.register('anaChannels', 'all', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Analysis channels to run.")
options.register('energyScales', 'all', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Event energy scales to run.")
options.register('productionMode', 'hh', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Selections that should be used for the production.")
options.register('tupleOutput', 'eventTuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Event tuple file.")
options.register('runSVfit', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Run SVfit algorithm on the selected tau pair.")
options.register('runKinFit', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Run HHKinFit algorithm for on the selected tau pair and all possible jet combinations.")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JSON file with lumi mask.")
options.parseArguments()

sampleConfig = importlib.import_module('h-tautau.Production.sampleConfig')
isData = sampleConfig.IsData(options.sampleType)

processName = 'tupleProduction'
process = cms.Process(processName)
process.options = cms.untracked.PSet()
process.options.wantSummary = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = options.globalTag
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(), eventsToProcess = cms.untracked.VEventRange('273158:247:387485495',
                '273158:248:389655799','273158:249:391412893','273158:270:423033611','273158:280:439046580',
				'273158:281:440165171','273158:282:441858587','273158:366:566581502','273158:441:676817062',
				'273158:442:678132807','273158:498:759027053','273158:566:854528991','273158:566:854633735',
				'273158:568:857031093','273158:571:862060594','273158:574:866112592','273158:693:1027720130',
				'273158:694:1028305820','273158:755:1109405653','273158:761:1116600141','273158:784:1146117685',
				'273158:789:1153712701','273158:792:1156869236','273158:793:1158025597','273158:870:1258977370',
				'273158:872:1261128991','273158:876:1266405514','273158:878:1269603134','273158:879:1270545938',
				'273158:879:1271053136','273158:880:1271905340','273158:885:1278228768','273158:929:1335163141',
				'273158:932:1338188812','273158:940:1349381100','273158:941:1350246859','273158:947:1357938046',
				'273158:951:1363097463','273158:958:1371274764','273158:1068:1510293276','273158:1075:1519530854','273158:1087:1534094776',
				'273158:1133:1590982760','273158:1138:1598122864','273158:1139:1598811935','273158:1139:1599132130'
                ))
process.TFileService = cms.Service('TFileService', fileName = cms.string(options.tupleOutput) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(0) )

if len(options.fileList) > 0:
    from AnalysisTools.Run.readFileList import *
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
    process.maxEvents.input = options.maxEvents

if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

process.load('RecoMET.METProducers.METSignificance_cfi')
process.load('RecoMET.METProducers.METSignificanceParams_cfi')

### re-run JEC and and use corrected jets for MET significance, if requested
if options.ReRunJEC:
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag('slimmedJets'),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
        payload = 'AK4PFchs' # Make sure to choose the appropriate levels and payload here!
    )
    if isData:
      process.patJetCorrFactorsReapplyJEC.levels.append('L2L3Residual')

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsReapplyJEC'))
    )

    process.METSignificance.srcPfJets = cms.InputTag('patJetsReapplyJEC', '', processName)
    process.METSignificance.srcMet = cms.InputTag('patpfMETT1', '', processName)

    process.JECsequence = cms.Sequence(process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC)
    JetCollectionName = 'patJetsReapplyJEC'
    #JetCollectionName = cms.InputTag('patJetsReapplyJEC', '', processName)
else:
   process.JECsequence = cms.Sequence()
   JetCollectionName = 'slimmedJets'
   #JetCollectionTag = cms.InputTag('slimmedJets')

from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
runMVAMET( process, jetCollectionPF = JetCollectionName )
process.MVAMET.srcLeptons  = cms.VInputTag('slimmedMuons', 'slimmedElectrons', 'slimmedTaus')
process.MVAMET.requireOS = cms.bool(False)

## Load module for Electron MVA ID
## It will append a value maps the miniAOD, that it's accesible throught a well Handle
## Example code here:
##  https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDwithMVADemo.cc#L99
## process.load('RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi')
##-------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD) ##also compute a maps with the electrons that pass an MVA cut
switchOnVIDElectronIdProducer(process, DataFormat.AOD)

# define which IDs we want to produce
id_modules = [ 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
               'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
               'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff' ]

#add them to the VID producer
for idmod in id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
#------------

process.tupleProductionSequence = cms.Sequence()
if options.anaChannels == 'all':
    channels = [ 'eTau', 'muTau', 'tauTau' ]
else:
    channels = re.split(',', options.anaChannels)

if options.energyScales == 'all':
    energyScales = [ 'Central', 'TauUp', 'TauDown', 'JetUp', 'JetDown' ]
else:
    energyScales = re.split(',', options.energyScales)

for channel in channels:
    producerName = 'tupleProducer_{}'.format(channel)
    producerClassName = 'TupleProducer_{}'.format(channel)
    hltPaths = sampleConfig.GetHltPaths(channel, options.sampleType)
    setattr(process, producerName, cms.EDAnalyzer(producerClassName,
        genParticles            = cms.InputTag('genParticles'),
        electronSrc             = cms.InputTag('slimmedElectrons'),
        eleTightIdMap           = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
        eleMediumIdMap          = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
        eleCutBasedVeto         = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
        tauSrc                  = cms.InputTag('slimmedTaus'),
        muonSrc                 = cms.InputTag('slimmedMuons'),
        vtxSrc                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
        jetSrc                  = cms.InputTag(JetCollectionName),
        fatJetSrc               = cms.InputTag('slimmedJetsAK8'),
        PUInfo                  = cms.InputTag('slimmedAddPileupInfo'),
        pfMETSrc                = cms.InputTag('slimmedMETs'),
        bits                    = cms.InputTag('TriggerResults', '', 'HLT'),
        prescales               = cms.InputTag('patTrigger'),
        objects                 = cms.InputTag('selectedPatTrigger'),
        metCov                  = cms.InputTag('METSignificance', 'METCovariance'),
        lheEventProducts        = cms.InputTag('externalLHEProducer'),
        genEventInfoProduct     = cms.InputTag('generator'),
        pruned                  = cms.InputTag('prunedGenParticles'),
        l1JetParticleProduct    = cms.InputTag('l1extraParticles', 'IsoTau'),
        isMC                    = cms.bool(not isData),
        applyTriggerMatch       = cms.bool(options.applyTriggerMatch),
        hltPaths                = cms.vstring(hltPaths),
        runSVfit                = cms.bool(options.runSVfit),
        runKinFit               = cms.bool(options.runKinFit),
        energyScales            = cms.vstring(energyScales),
        productionMode          = cms.string(options.productionMode)
    ))
    process.tupleProductionSequence += getattr(process, producerName)

process.p = cms.Path(
    process.egmGsfElectronIDSequence *
    process.electronMVAValueMapProducer *
    process.JECsequence *
    process.METSignificance *
    process.tupleProductionSequence
)

#print process.dumpPython()
