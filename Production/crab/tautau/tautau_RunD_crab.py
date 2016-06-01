from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'tautau_RunD'
config.General.workArea = '2015Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/Production_76X_tautau.py'
config.JobType.pyCfgParams = ['sampleType=Run2015D','isData=True','ReRunJEC=True','computeHT=False']

config.Data.inputDataset = '/Tau/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 3000
config.Data.lumiMask = '/afs/cern.ch/work/c/ccaputo/private/HHbbTauTau/Run2/DEV/CMSSW_7_6_3_patch2/src/h-tautau/Production/json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/76X_Production_v3/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
