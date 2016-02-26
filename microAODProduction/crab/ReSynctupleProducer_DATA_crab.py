from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'ReSyncTuple_Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=True','runOnCrab=True','sampleType=Run2015D']

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
    ##      RunD     ##
    ###################
    config.General.requestName = 'RunD'
    config.Data.inputDataset = '/SingleMuon/ccaputo-crab_RunD_3rd-c953508963d94912a54c4204017e3a7a/USER'
    config.Data.unitsPerJob = 8000
    submit(config)

    config.General.requestName = 'RunD_Oct'
    config.Data.inputDataset = '/SingleMuon/ccaputo-crab_RunD_Oct_3rd-c953508963d94912a54c4204017e3a7a/USER'
    config.Data.unitsPerJob = 8000
    submit(config)
