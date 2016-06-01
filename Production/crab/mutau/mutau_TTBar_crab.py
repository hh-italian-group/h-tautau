from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'mutau_TTbar_ext3'
config.General.workArea = 'TTbar'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/Production_76X_mutau.py'
config.JobType.pyCfgParams = ['sampleType=Fall15MC','isData=False','ReRunJEC=True','computeHT=False']

config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/76X_Production_v3/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
