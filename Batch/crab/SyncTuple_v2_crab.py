# Produce SyncTuple v2.
# This file is part of https://github.com/hh-italian-group/h-tautau.

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SUSYGluGluToHToTauTau_M160_PAT_Sync'
config.General.workArea = 'SUSYGluGluToHToTauTau_M160'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'

config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 2000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/160108_PAT_Sync/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
