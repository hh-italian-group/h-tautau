# Produce Run D data.
# This file is part of https://github.com/hh-italian-group/h-tautau.

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'RunD_SyncTree_v2'
config.General.workArea = '2015Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=True','runOnCrab=True','sampleType=Run2015D']

config.Data.inputDataset = '/SingleMuon/ccaputo-crab_RunD_micro_Sync_v2-c953508963d94912a54c4204017e3a7a/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 8000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/RunOnPublishedDataset/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
