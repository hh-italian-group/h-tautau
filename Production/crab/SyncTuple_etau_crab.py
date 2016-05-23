from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ETau_tauES_GenMatch'
config.General.workArea = 'SUSYGluGluToHToTauTau_M160'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/Production_76X_etau.py'
config.JobType.pyCfgParams = ['isData=False','ReRunJEC=True','sampleType=Fall15MC','computeHT=False']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['microAOD.root','syncTree.root','etau_cuts.root']
#config.JobType.outputFiles = ['microAOD.root','syncTree.root','mutau_cuts.root']


config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 4000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/Sync_160522/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
