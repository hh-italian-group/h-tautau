from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'ReSyncTuple_MC'
config.General.requestName = 'GluGluToRadionToHHTo2B2Tau_M-250'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=False','runOnCrab=True','computeHT=false']

config.Data.inputDataset = '/GluGluToRadionToHHTo2B2Tau_M-250_narrow_13TeV-madgraph/ccaputo-crab_GluGluToRadionToHHTo2B2Tau_M-250-46ffc704b6280aca3e8f5a327bcc5603/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/ReSyncTuple/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'

#if __name__ == '__main__':

#    from CRABAPI.RawCommand import crabCommand
#    from CRABClient.ClientExceptions import ClientException
#    from httplib import HTTPException


#    def submit(config):
#        try:
#            crabCommand('submit', config = config)
#        except HTTPException as hte:
#            print "Failed submitting task: %s" % (hte.headers)
#        except ClientException as cle:
#            print "Failed submitting task: %s" % (cle)

#    ###################
#    ##      Hhh      ##
#    ###################
