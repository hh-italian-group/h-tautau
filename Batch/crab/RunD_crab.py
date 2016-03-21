# Produce Run D data.
# This file is part of https://github.com/hh-italian-group/h-tautau.

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'RunD'
config.General.workArea = '2015Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'
config.JobType.pyCfgParams = ['sampleType=Run2015D','isData=True','computeHT=False']

config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 10000
config.Data.lumiMask = '/cmshome/caputo/HH_bbTauTau/Run2/CMSSW_7_6_3_patch2/src/HHbbTauTau/microAODProduction/json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/76X_Production_v1/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T3_IT_MIB'
