from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'GluGluToRadionToHHTo2B2Tau'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'
config.JobType.pyCfgParams = ['sampleType=Fall15MC','isData=False','computeHT=False']

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/76X_Production_v1/' # or '/store/group/<subdir>'
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

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-250'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-250_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-260'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-260_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-270'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-270_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-280'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-280_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-300'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-320'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-320_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-340'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-340_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-350'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-350_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-400'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-400_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-450'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-450_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-500'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-500_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-550'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-550_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-600'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-600_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-650'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-650_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-700'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-700_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-800'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-800_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)

    config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-900'
    config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-900_narrow_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1000
    submit(config)



