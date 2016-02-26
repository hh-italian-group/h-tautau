from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'ReSyncTuple_MC'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=False','runOnCrab=True']

config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/ReSyncTuple/' # or '/store/group/<subdir>'
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
    ##    DiBoson    ##
    ###################
    config.General.requestName = 'VVTo2L2Nu_13TeV'
    config.Data.inputDataset = '/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_VVTo2L2Nu_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 3000
    submit(config)

    config.General.requestName = 'ZZTo2L2Q_13TeV'
    config.Data.inputDataset = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_ZZTo2L2Q_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 5000
    submit(config)

    config.General.requestName = 'ZZTo4L_13TeV'
    config.Data.inputDataset = '/ZZTo4L_13TeV-amcatnloFXFX-pythia8/ccaputo-crab_ZZTo4L_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 5000
    submit(config)

    config.General.requestName = 'WWTo1L1Nu2Q_13TeV'
    config.Data.inputDataset = '/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_WWTo1L1Nu2Q_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 2000
    submit(config)

    config.General.requestName = 'WZTo2L2Q_13TeV'
    config.Data.inputDataset = '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_WZTo2L2Q_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 9000
    submit(config)

    config.General.requestName = 'WZJets_13TeV'
    config.Data.inputDataset = '/WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ccaputo-crab_WZJets_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 5000
    submit(config)

    config.General.requestName = 'WZTo1L3Nu_13TeV'
    config.Data.inputDataset = '/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_WZTo1L3Nu_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 700
    submit(config)

    config.General.requestName = 'WZTo1L1Nu2Q_13TeV'
    config.Data.inputDataset = '/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/ccaputo-crab_WZTo1L1Nu2Q_13TeV-2470bbdf15448d116b7dc09c605c0dbb/USER'
    config.Data.unitsPerJob = 7000
    submit(config)

    ###################
    ##     TTBar     ##
    ###################
    config.General.requestName = 'TTbar_ext3'
    config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ccaputo-crab_TTbar_ext3_3rd-46ffc704b6280aca3e8f5a327bcc5603/USER'
    config.Data.unitsPerJob = 40000
    submit(config)

#    ###################
#    ##     WJets     ##
#    ###################
#    config.General.requestName = 'WJetsToLNu'
#    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ccaputo-crab_WJetsToLNu_3rd-46ffc704b6280aca3e8f5a327bcc5603/USER'
#    config.Data.unitsPerJob = 50000
#    submit(config)

#    ###################
#    ##     DYJets    ##
#    ###################
#    config.General.requestName = 'DYJetsToLL_M-50'
#    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ccaputo-crab_DYJetsToLL_M-50_3rd-46ffc704b6280aca3e8f5a327bcc5603/USER'
#    config.Data.unitsPerJob = 15000
#    submit(config)
