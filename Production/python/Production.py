# Produce EventTuple for all channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Global Tag to use.")
options.register('sampleType', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
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
options.register('applyRecoilCorr', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply Met Recoil Corrections")
options.register('nJetsRecoilCorr', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                        "Number of Additional Jets for Recoil Correction")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "JSON file with lumi mask.")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "List of events to process.")
options.register('saveGenTopInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for top quarks.")
options.register('saveGenBosonInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for bosons.")
options.register('saveGenJetInfo', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for jets.")
options.register('dumpPython', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Dump full config into stdout.")
options.register('numberOfThreads', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                        "Number of threads.")

options.parseArguments()

sampleConfig = importlib.import_module('h-tautau.Production.sampleConfig')
isData = sampleConfig.IsData(options.sampleType)
period = sampleConfig.GetPeriod(options.sampleType)
triggerCfg = sampleConfig.GetTriggerCfg(period)

processName = 'tupleProduction'
process = cms.Process(processName)
process.options = cms.untracked.PSet()
process.options.wantSummary = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)
process.options.numberOfThreads = cms.untracked.uint32(options.numberOfThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.GlobalTag.globaltag = options.globalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.TFileService = cms.Service('TFileService', fileName = cms.string(options.tupleOutput) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(0) )

if len(options.fileList) > 0:
    from AnalysisTools.Run.readFileList import *
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
    process.maxEvents.input = options.maxEvents

if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

## and add btag discriminator to the event content
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

if period == 'Run2016':
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
        btagDiscriminators = [
            'pfDeepCSVJetTags:probudsg',
            'pfDeepCSVJetTags:probb',
            'pfDeepCSVJetTags:probc',
            'pfDeepCSVJetTags:probbb',
        ]
    )

if period == 'Run2017':
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
    )
    process.jecSequence = cms.Sequence(process.patJetCorrFactors * process.updatedPatJets * process.selectedUpdatedPatJets)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData = isData)

MetInputTag = cms.InputTag('slimmedMETs', '', processName)

### MET filters for 2016
if period == 'Run2016':
    process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
    process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
    process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

    process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
    process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
    process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

## Load module for Electron MVA ID
## It will append a value maps the miniAOD, that it's accesible throught a well Handle
## Example code here:
##  https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDwithMVADemo.cc#L99
process.load('RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi')
##-------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD) ##also compute a maps with the electrons that pass an MVA cut
#switchOnVIDElectronIdProducer(process, DataFormat.AOD)

# define which IDs we want to produce

if period == 'Run2016':
    id_modules = [ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff' ]


if period == 'Run2017':
    id_modules = [ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                   'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff' ]

#add them to the VID producer
for idmod in id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

if period == 'Run2016':
    tauSrc_InputTag = cms.InputTag('slimmedTaus')

if period == 'Run2017':
    import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
    updatedTauName = "slimmedTausNewID"
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(
        process, cms, debug = True, updatedTauName = updatedTauName,
        toKeep = [ "2017v2", "dR0p32017v2", "newDM2017v2", "2016v1", "newDM2016v1",
                   "deepTau2017v1", "DPFTau_2016_v0", "againstEle2018", ]
    )
    tauIdEmbedder.runTauID()
    tauSrc_InputTag = cms.InputTag('slimmedTausNewID')

if period == 'Run2016':
    tauAntiEle = importlib.import_module('h-tautau.Production.runTauAgainstElectron')
    tauAntiEle.rerunAgainstElectron(process, process.NewTauIDsEmbedded)


### Top gen level info
process.topGenSequence = cms.Sequence()
if options.saveGenTopInfo:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    process.decaySubset.fillMode = cms.string("kME")
    process.initSubset.src = cms.InputTag('prunedGenParticles')
    process.decaySubset.src = cms.InputTag('prunedGenParticles')
    process.decaySubset.runMode = cms.string("Run2")
    process.topGenSequence += process.makeGenEvt

if options.anaChannels == 'all':
    channels = [ 'eTau', 'muTau', 'tauTau' ]
else:
    channels = re.split(',', options.anaChannels)

if options.energyScales == 'all':
    energyScales = [ 'Central', 'TauUp', 'TauDown', 'JetUp', 'JetDown' ]
else:
    energyScales = re.split(',', options.energyScales)

### Tuple production sequence

if period == 'Run2016':
    eleTightIdMap_InputTag           = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80')
    eleMediumIdMap_InputTag          = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90')
#    tauSrc_InputTag                  = cms.InputTag('slimmedTaus')
    objects_InputTag                 = cms.InputTag('selectedPatTrigger')

if period == 'Run2017':
    eleTightIdMap_InputTag           = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80')
    eleMediumIdMap_InputTag          = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90')
#    tauSrc_InputTag                  = cms.InputTag('NewTauIDsEmbedded')
    objects_InputTag                 = cms.InputTag('slimmedPatTrigger')

process.summaryTupleProducer = cms.EDAnalyzer('SummaryProducer',
    isMC            = cms.bool(not isData),
    saveGenTopInfo  = cms.bool(options.saveGenTopInfo),
    lheEventProduct = cms.InputTag('externalLHEProducer'),
    genEvent        = cms.InputTag('generator'),
    topGenEvent     = cms.InputTag('genEvt'),
    puInfo          = cms.InputTag('slimmedAddPileupInfo'),
    triggerCfg      = cms.string(triggerCfg),
    channels        = cms.vstring(channels)
)

process.tupleProductionSequence = cms.Sequence(process.summaryTupleProducer)


for channel in channels:
    producerName = 'tupleProducer_{}'.format(channel)
    producerClassName = 'TupleProducer_{}'.format(channel)

    setattr(process, producerName, cms.EDAnalyzer(producerClassName,
        electronSrc             = cms.InputTag('slimmedElectrons'),
        eleTightIdMap           = eleTightIdMap_InputTag,
        eleMediumIdMap          = eleMediumIdMap_InputTag,
        tauSrc                  = tauSrc_InputTag,
        muonSrc                 = cms.InputTag('slimmedMuons'),
        vtxSrc                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
        jetSrc                  = cms.InputTag('selectedUpdatedPatJets'),
        fatJetSrc               = cms.InputTag('slimmedJetsAK8'),
        PUInfo                  = cms.InputTag('slimmedAddPileupInfo'),
        pfMETSrc                = MetInputTag,
        prescales               = cms.InputTag('patTrigger'),
        objects                 = objects_InputTag,
        lheEventProducts        = cms.InputTag('externalLHEProducer'),
        genEventInfoProduct     = cms.InputTag('generator'),
        topGenEvent             = cms.InputTag('genEvt'),
        genParticles            = cms.InputTag('prunedGenParticles'),
        genJets                 = cms.InputTag('slimmedGenJets'),
        l1JetParticleProduct    = cms.InputTag('l1extraParticles', 'IsoTau'),
        isMC                    = cms.bool(not isData),
        applyTriggerMatch       = cms.bool(options.applyTriggerMatch),
        runSVfit                = cms.bool(options.runSVfit),
        runKinFit               = cms.bool(options.runKinFit),
        applyRecoilCorr         = cms.bool(options.applyRecoilCorr),
        nJetsRecoilCorr         = cms.int32(options.nJetsRecoilCorr),
        energyScales            = cms.vstring(energyScales),
        productionMode          = cms.string(options.productionMode),
        period                  = cms.string(period),
        triggerCfg              = cms.string(triggerCfg),
        saveGenTopInfo          = cms.bool(options.saveGenTopInfo),
        saveGenBosonInfo        = cms.bool(options.saveGenBosonInfo),
        saveGenJetInfo          = cms.bool(options.saveGenJetInfo),
        rho                     = cms.InputTag('fixedGridRhoAll'),
    ))

    if period == 'Run2016':
        getattr(process, producerName).badPFMuonFilter = cms.InputTag('BadPFMuonFilter')
        getattr(process, producerName).badChCandidateFilter = cms.InputTag('BadChargedCandidateFilter')

    process.tupleProductionSequence += getattr(process, producerName)

if period == 'Run2016':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.electronMVAValueMapProducer *
        process.fullPatMetSequence *
        process.BadPFMuonFilter *
        process.BadChargedCandidateFilter *
        process.topGenSequence *
        process.rerunMvaIsolationSequence *
        process.rerunDiscriminationAgainstElectronMVA6 *
        process.tupleProductionSequence
    )

if period == 'Run2017':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.electronMVAValueMapProducer *
        process.jecSequence *
        process.rerunMvaIsolationSequence *
        getattr(process, updatedTauName) *
        process.fullPatMetSequence *
        process.topGenSequence *
        process.tupleProductionSequence
    )

if options.dumpPython:
    print process.dumpPython()
