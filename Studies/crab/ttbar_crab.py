from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ttBar_Rejection_GENStudies_FirstProduction'
config.General.workArea = 'ttBar_Rejection_GENStudies'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/ttbarGenStudies_cfg.py'
# config.JobType.pyCfgParams = ['isData=False','ReRunJEC=True','sampleType=Fall15MC','computeHT=False']
# config.JobType.disableAutomaticOutputCollection = True
# config.JobType.outputFiles = ['microAOD.root','syncTree.root','mutau_cuts.root','etau_cuts.root','tautau_cuts.root']
# #config.JobType.outputFiles = ['microAOD.root','syncTree.root','mutau_cuts.root']


config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 30000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
