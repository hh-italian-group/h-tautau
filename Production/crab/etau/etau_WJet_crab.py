from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'etau_WJetsToLNu_HT'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/Production_76X_etau.py'
config.JobType.pyCfgParams = ['sampleType=Fall15MC','isData=False','ReRunJEC=True','computeHT=True']

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/76X_Production_v3/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    ###################
    ##     WJets     ##
    ###################
    config.General.requestName = 'WJetsToLNu'
    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 37000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-100to200'
    config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 8000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-200to400'
    config.Data.inputDataset = '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 7000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-400to600'
    config.Data.inputDataset = '/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 7000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-600toInf'
    config.Data.inputDataset = '/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 7000
    submit(config)
