from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'GluGluToRadionToHHTo2B2Tau'
config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-550'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_EleID.py'

config.Data.inputDBS = 'global'
config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-550_narrow_13TeV-madgraph/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.unitsPerJob = 1000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/FirstProduction/' # or '/store/group/<subdir>'
config.Data.publication = True

config.Site.storageSite = 'T2_IT_Bari'

