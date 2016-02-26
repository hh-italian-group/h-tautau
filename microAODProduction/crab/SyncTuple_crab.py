from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SUSYGluGluToHToTauTau_M160_76X_Pisa'
config.General.workArea = 'SUSYGluGluToHToTauTau_M160'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'
config.JobType.pyCfgParams = ['isData=False','sampleType=Fall15MC','computeHT=False']

config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 5000
#config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/160215_Sync_76x_v1/' # or '/store/group/<subdir>'
config.Data.outLFNDirBase = '/store/user/ccaputo/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Pisa'
