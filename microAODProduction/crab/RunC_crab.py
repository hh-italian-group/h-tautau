from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'RunC_micro_Sync_v2'
config.General.workArea = '2015Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_EleID_DATA.py'

config.Data.inputDataset = '/SingleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.lumiMask = '/cmshome/caputo/HH_bbTauTau/Run2/CMSSW_7_4_12_patch4/src/HHbbTauTau/microAODProduction/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 10000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/SecondProduction/' # or '/store/group/<subdir>'
config.Data.publication = True

config.Site.storageSite = 'T2_IT_Bari'
