# Produce EventTuple for all channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('sampleType', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
options.register('applyTriggerMatch', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger matching for signal objects. Default: True")
options.register('applyTriggerMatchCut', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger matching Cut not for signal samples. Default: True")
options.register('fileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "List of root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Prefix to add to input file names.")
options.register('anaChannels', 'all', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Analysis channels to run.")
options.register('tupleOutput', 'eventTuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Event tuple file.")
options.register('applyTriggerCut', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger cut for signal objects. Default: True")
options.register('storeLHEinfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Store LHE information. Default: False")
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
options.register('saveGenParticleInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for particles.")
options.register('isEmbedded', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Is DY embedded sample.")
options.register('dumpPython', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Dump full config into stdout.")
options.register('numberOfThreads', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                        "Number of threads.")

options.parseArguments()

sampleConfig = importlib.import_module('h-tautau.Production.sampleConfig')
isData = sampleConfig.IsData(options.sampleType)
period = sampleConfig.GetPeriod(options.sampleType)
triggerCfg = sampleConfig.GetTriggerCfg(options.sampleType)

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

process.GlobalTag.globaltag = sampleConfig.GetGlobalTag(options.sampleType)
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.TFileService = cms.Service('TFileService', fileName = cms.string(options.tupleOutput) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from AnalysisTools.Run.readFileList import *
if len(options.fileList) > 0:
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
elif len(options.inputFiles) > 0:
    addFilesToList(process.source.fileNames, options.inputFiles, options.fileNamePrefix)
if options.maxEvents > 0:
    process.maxEvents.input = options.maxEvents

if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

## and add btag discriminator to the event content
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools
btagVector = []

if period == 'Run2018':
    btagVector.append('None')

if period == 'Run2017':
    btagVector2017 = [
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probc',
        'pfDeepFlavourJetTags:probuds',
        'pfDeepFlavourJetTags:probg'
    ]
    btagVector.extend(btagVector2017)

if period == 'Run2016':
    btagVector2016 = [
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probc',
        'pfDeepFlavourJetTags:probuds',
        'pfDeepFlavourJetTags:probg',
        'pfDeepCSVJetTags:probudsg',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb'
    ]
    btagVector.extend(btagVector2016)


jec_levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring(jec_levels), 'None'),
    btagDiscriminators = btagVector,
    postfix='NewDFTraining'
)
if period == 'Run2016':
    process.jecSequence = cms.Sequence(process.patJetCorrFactorsNewDFTraining *
                                       process.updatedPatJetsNewDFTraining *
                                       process.patJetCorrFactorsTransientCorrectedNewDFTraining *
                                       process.pfImpactParameterTagInfosNewDFTraining *
                                       process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining *
                                       process.pfDeepCSVTagInfosNewDFTraining *
                                       process.pfDeepCSVJetTagsNewDFTraining *
                                       process.pfDeepFlavourTagInfosNewDFTraining *
                                       process.pfDeepFlavourJetTagsNewDFTraining *
                                       process.updatedPatJetsTransientCorrectedNewDFTraining *
                                       process.selectedUpdatedPatJetsNewDFTraining)

if period == 'Run2017':
    process.jecSequence = cms.Sequence(process.patJetCorrFactorsNewDFTraining *
                                       process.updatedPatJetsNewDFTraining *
                                       process.patJetCorrFactorsTransientCorrectedNewDFTraining *
                                       process.pfImpactParameterTagInfosNewDFTraining *
                                       process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining *
                                       process.pfDeepCSVTagInfosNewDFTraining *
                                       process.pfDeepFlavourTagInfosNewDFTraining *
                                       process.pfDeepFlavourJetTagsNewDFTraining *
                                       process.updatedPatJetsTransientCorrectedNewDFTraining *
                                       process.selectedUpdatedPatJetsNewDFTraining)

if period == 'Run2018':
    process.jecSequence = cms.Sequence(process.patJetCorrFactorsNewDFTraining *
                                       process.updatedPatJetsNewDFTraining *
                                       #process.patJetCorrFactorsTransientCorrectedNewDFTraining *
                                       #process.pfImpactParameterTagInfosNewDFTraining *
                                       #process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining *
                                       #process.updatedPatJetsTransientCorrectedNewDFTraining *
                                       process.selectedUpdatedPatJetsNewDFTraining)


### MET filters for 2016 and MET recipe for 2016
if period == 'Run2016':
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process, isData = isData or options.isEmbedded)
    MetInputTag = cms.InputTag('slimmedMETs', '', processName)

customMetFilters = cms.PSet()
### MET filters for 2017 and MET recipe for 2017
if period == 'Run2017' or period == 'Run2018':
    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

    baddetEcallist = cms.vuint32(
        [872439604,872422825,872420274,872423218,
        872423215,872416066,872435036,872439336,
        872420273,872436907,872420147,872439731,
        872436657,872420397,872439732,872439339,
        872439603,872422436,872439861,872437051,
        872437052,872420649,872421950,872437185,
        872422564,872421566,872421695,872421955,
        872421567,872437184,872421951,872421694,
        872437056,872437057,872437313,872438182,
        872438951,872439990,872439864,872439609,
        872437181,872437182,872437053,872436794,
        872436667,872436536,872421541,872421413,
        872421414,872421031,872423083,872421439])


    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
        "EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt        = cms.double(50.),
        baddetEcal    = baddetEcallist,
        taggingMode = cms.bool(True),
        debug = cms.bool(False)
        )
    customMetFilters.ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")

if period == 'Run2017':
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

    runMetCorAndUncFromMiniAOD (
            process,
            isData = isData or options.isEmbedded, # false for MC
            fixEE2017 = True,
            fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
            postfix = "ModifiedMET"
    )

    MetInputTag = cms.InputTag('slimmedMETsModifiedMET', '', processName)

if period == 'Run2018':
    MetInputTag = cms.InputTag('slimmedMETs')


# Reweighting recipe to emulate Level 1 ECAL prefiring from https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
if period == 'Run2016' or period == 'Run2017':
    from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
    prefiring_eras = { 'Run2016': "2016BtoH", 'Run2017': "2017BtoF" }
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(DataEra = cms.string(prefiring_eras[period]),
                                                                  UseJetEMPt = cms.bool(False),
                                                                  PrefiringRateSystematicUncty = cms.double(0.2),
                                                                  SkipWarnings = False)

# Update electron ID following recommendations from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
ele_era = { 'Run2016': '2016-Legacy', 'Run2017': '2017-Nov17ReReco', 'Run2018': '2018-Prompt'}
setupEgammaPostRecoSeq(process, runVID=True, runEnergyCorrections=False, era=ele_era[period])

import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
updatedTauName = "slimmedTausNewID"
tauIdEmbedder = tauIdConfig.TauIDEmbedder(
    process, cms, debug = True, updatedTauName = updatedTauName,
    toKeep = [ "2017v2", "deepTau2017v2p1",  "againstEle2018", "2016v1"]
)
tauIdEmbedder.runTauID()
tauSrc_InputTag = cms.InputTag('slimmedTausNewID')

jetSrc_InputTag                  = cms.InputTag('selectedUpdatedPatJetsNewDFTraining')
objects_InputTag                 = cms.InputTag('slimmedPatTrigger')

# Update pileup jet ID following recommendations from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId, _chsalgos_81x, _chsalgos_94x, _chsalgos_102x
pu_jet_id_wp = { 'Run2016': _chsalgos_81x, 'Run2017': _chsalgos_94x, 'Run2018': _chsalgos_102x}
process.updatedPileupJetId = pileupJetId.clone(
    jets = jetSrc_InputTag,
    inputIsCorrected = True,
    applyJec = False,
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    algos = cms.VPSet(pu_jet_id_wp[period]),
)

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
    channels = [ 'eTau', 'muTau', 'tauTau', 'muMu' ]
else:
    channels = re.split(',', options.anaChannels)

### Tuple production sequence

if period == 'Run2016' and options.isEmbedded :
    objects_InputTag             = cms.InputTag('selectedPatTrigger')

genJets_inputTag = cms.InputTag('slimmedGenJets')

if period == 'Run2017' and options.isEmbedded :
    genJets_inputTag = cms.InputTag('slimmedGenJetsAK8SoftDropSubJets')

if period == 'Run2018' and options.isEmbedded :
    genJets_inputTag = cms.InputTag('slimmedGenJetsAK8SoftDropSubJets')


process.summaryTupleProducer = cms.EDAnalyzer('SummaryProducer',
    isMC            = cms.bool(not isData or options.isEmbedded),
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
        tauSrc                  = tauSrc_InputTag,
        muonSrc                 = cms.InputTag('slimmedMuons'),
        vtxSrc                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
        jetSrc                  = jetSrc_InputTag,
        fatJetSrc               = cms.InputTag('slimmedJetsAK8'),
        PUInfo                  = cms.InputTag('slimmedAddPileupInfo'),
        pfMETSrc                = MetInputTag,
        genMetSrc               = cms.InputTag('genMetTrue'),
        prescales               = cms.InputTag('patTrigger'),
        objects                 = objects_InputTag,
        lheEventProducts        = cms.InputTag('externalLHEProducer'),
        genEventInfoProduct     = cms.InputTag('generator'),
        topGenEvent             = cms.InputTag('genEvt'),
        genParticles            = cms.InputTag('prunedGenParticles'),
        genJets                 = genJets_inputTag,
        l1JetParticleProduct    = cms.InputTag('l1extraParticles', 'IsoTau'),
        isMC                    = cms.bool(not isData or options.isEmbedded),
        applyTriggerMatch       = cms.bool(options.applyTriggerMatch),
        applyTriggerMatchCut    = cms.bool(options.applyTriggerMatchCut),
        applyTriggerCut         = cms.bool(options.applyTriggerCut),
        storeLHEinfo            = cms.bool(options.storeLHEinfo),
        nJetsRecoilCorr         = cms.int32(options.nJetsRecoilCorr),
        period                  = cms.string(period),
        triggerCfg              = cms.string(triggerCfg),
        saveGenTopInfo          = cms.bool(options.saveGenTopInfo),
        saveGenBosonInfo        = cms.bool(options.saveGenBosonInfo),
        saveGenJetInfo          = cms.bool(options.saveGenJetInfo),
        saveGenParticleInfo     = cms.bool(options.saveGenParticleInfo),
        isEmbedded              = cms.bool(options.isEmbedded),
        rho                     = cms.InputTag('fixedGridRhoAll'),
        customMetFilters        = customMetFilters,
        updatedPileupJetIdDiscr = cms.InputTag('updatedPileupJetId', 'fullDiscriminant'),
        updatedPileupJetId      = cms.InputTag('updatedPileupJetId', 'fullId'),
    ))

    process.tupleProductionSequence += getattr(process, producerName)

if period == 'Run2016':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.egammaPostRecoSeq *
        process.jecSequence *
        process.updatedPileupJetId *
        process.rerunMvaIsolationSequence *
        getattr(process, updatedTauName) *
        process.fullPatMetSequence *
        process.topGenSequence *
        process.prefiringweight *
        process.tupleProductionSequence
    )

if period == 'Run2017':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.egammaPostRecoSeq *
        process.jecSequence *
        process.updatedPileupJetId *
        process.rerunMvaIsolationSequence *
        getattr(process, updatedTauName) *
        process.fullPatMetSequenceModifiedMET *
        process.ecalBadCalibReducedMINIAODFilter *
        process.topGenSequence *
        process.prefiringweight*
        process.tupleProductionSequence
    )

if period == 'Run2018':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.egammaPostRecoSeq *
        process.jecSequence *
        process.updatedPileupJetId *
        process.rerunMvaIsolationSequence *
        getattr(process, updatedTauName) *
        process.ecalBadCalibReducedMINIAODFilter *
        process.topGenSequence *
        process.tupleProductionSequence
    )

if options.dumpPython:
    print process.dumpPython()
