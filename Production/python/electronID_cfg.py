# Electron ID configuration.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root',
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00610EE7-C213-E511-842C-00304833529A.root',
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00623DCC-A813-E511-A302-0025905B85A2.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root',
    )

personalFileMiniAOD = cms.untracked.vstring('file:768F5AFB-D771-E511-9ABD-B499BAABD280.root') ##skimBranches_Signal_miniV2.root
DYSample = cms.untracked.vstring("/store/user/ccaputo/HHbbtautau/Run2/DYSample_forHT.root")
RunD_mini = cms.untracked.vstring('file:/lustre/cms//store/user/ccaputo/HHbbtautau/Run2/FirstProduction/SingleMuon/crab_RunD/160104_205440/0000/skimBranches_Data_140.root')
# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = DYSample
    outputFile = "tauID_ntuple.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles,
#                              eventsToProcess = cms.untracked.VEventRange('1:97841')
#                             eventsToProcess = cms.untracked.VEventRange(
#                                            '1:4854','1:5650','1:12154', '1:15313', '1:18377', '1:22017', '1:22477',
#                                            '1:30133', '1:30975', '1:32636', '1:34615', '1:48908', '1:50869', '1:50980',
#                                            '1:53592', '1:58895', '1:59329', '1:61403', '1:62118', '1:64539', '1:65259',
#                                            '1:68675', '1:69486', '1:74288', '1:76282', '1:77713', '1:78710', '1:82726',
#                                            '1:85391', '1:89033', '1:94331', '1:96211', '1:96265', '1:97841', '1:98466',
#                                            '1:101778','1:106256','1:108005','1:115654','1:118548','1:124470','1:125911',
#                                            '1:128080','1:133973','1:145636','1:145782','1:148435','1:153388','1:154989',
#                                            '1:155222','1:155606','1:157896','1:158203','1:159524','1:159677','1:164960',
#                                            '1:173976','1:176806','1:178151','1:178506','1:185060','1:185621','1:206157',
#                                            '1:218844','1:224633','1:225924','1:226354','1:228458','1:229620','1:233295',
#                                            '1:235780','1:236624','1:240021','1:248476','1:248626','1:252329','1:261874',
#                                            '1:269397','1:271785','1:273917','1:275016','1:276244','1:279778','1:282814',
#                                            '1:285986','1:287025','1:296181','1:306248','1:311128','1:312811','1:313864',
#                                            '1:319594','1:321652','1:337202','1:354113','1:354826','1:356172','1:358279',
#                                            '1:361468','1:368865','1:372553','1:374591','1:374834','1:376422','1:377092',
#                                            '1:377857','1:381228','1:382663','1:389025','1:389619','1:395620','1:396061',
#                                            '1:403741','1:404549','1:407695','1:411717','1:411885','1:415110','1:421456',
#                                            '1:439435','1:440533','1:448589','1:449538','1:452435','1:458941','1:461904',
#                                            '1:475696','1:476634','1:481900','1:485343','1:485435','1:488045','1:488544',
#                                            '1:489345','1:489596','1:489812','1:494968','1:495099')
                                            )

#
# Set up electron ID (VID framework)
#

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## turn on VID producer, indicate data format  to be
## DataFormat.AOD or DataFormat.MiniAOD, as appropriate
#if useAOD == True :
#    dataFormat = DataFormat.AOD
#else :
#    dataFormat = DataFormat.MiniAOD

#switchOnVIDElectronIdProducer(process, dataFormat)

## define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

##add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#
# Configure an example module for user analysis with electrons
#
print "="*50
print "Setting the ntupler"
print "="*50

process.synctupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc    = cms.InputTag("slimmedTaus"),
                                 muonSrc   = cms.InputTag("slimmedMuons"),
                                 vtxSrc    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 pfMETSrc  = cms.InputTag("slimmedMETsNoHF"),
                                 jetSrc    = cms.InputTag("slimmedJets"),
                                 bits      = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects   = cms.InputTag("selectedPatTrigger"),
                                 sampleType = cms.string("Spring15MC"),
                                )

process.ntupler = cms.EDAnalyzer('ElectronsIDAnalyzer',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 # ... none ...
                                 #
                                 electronBool = cms.untracked.bool(False), ##Switch btw electron and tau ntupler
                                                                           ## TO BE IMPROVED
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 tauSrc  = cms.InputTag("slimmedTaus"),
                                 #
                                 # ID decisions (common to all formats)
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                 #
                                 # ValueMaps with MVA results
                                 #
                                 mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                 mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories")
                                )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

print "="*50
print "Starting the ntupler"
print "="*50

# Make sure to add the ID sequence upstream from the user analysis module
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
process.p = cms.Path(process.synctupler)
