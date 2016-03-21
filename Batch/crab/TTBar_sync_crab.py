# Produce ttbar sync.
# This file is part of https://github.com/hh-italian-group/h-tautau.

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TTbar_ext3_SyncTree'
config.General.workArea = 'TTbar_SyncTree'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=False','runOnCrab=True']

config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ccaputo-crab_TTbar_ext3-2470bbdf15448d116b7dc09c605c0dbb/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/RunOnPublishedDataset/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
