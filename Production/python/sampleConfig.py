# Configurations dependent on the sample type.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import sys
from sets import Set
import FWCore.ParameterSet.Config as cms

mcSampleTypes = Set([ 'MC_16', 'MC_17', 'MC_18', 'Emb_16', 'Emb_17', 'Emb_18ABC', 'Emb_18D' ])
dataSampleTypes = Set([ 'Run2016' , 'Run2017', 'Run2018ABC', 'Run2018D' ])

periodDict = { 'MC_16' : 'Run2016',
               'Run2016' : 'Run2016',
               'Emb_16' : 'Run2016',
               'MC_17' : 'Run2017',
               'Run2017' : 'Run2017',
               'Emb_17' : 'Run2017',
               'MC_18' : 'Run2018',
               'Run2018ABC' : 'Run2018',
               'Run2018D' : 'Run2018',
               'Emb_18ABC' : 'Run2018',
               'Emb_18D' : 'Run2018'
             }

globalTagMap = { 'MC_16' : '102X_mcRun2_asymptotic_v7',
                 'Run2016' : '102X_dataRun2_v12',
                 'Emb_16' : '102X_dataRun2_v12',
                 #'Emb_16' : '80X_dataRun2_2016SeptRepro_v7',
                 'MC_17' : '102X_mc2017_realistic_v7',
                 'Run2017' : '102X_dataRun2_v12',
                 'Emb_17' : '102X_dataRun2_v12',
                 'MC_18' : '102X_upgrade2018_realistic_v20',
                 'Run2018ABC' : '102X_dataRun2_v12',
                 'Run2018D' : '102X_dataRun2_Prompt_v15',
                 'Emb_18ABC' : '102X_dataRun2_v12',
                 'Emb_18D' : '102X_dataRun2_Prompt_v15'
               }

hltPaths = {
    'MC_16'   : 'h-tautau/Production/data/triggers_2016.cfg',
    'Run2016'   : 'h-tautau/Production/data/triggers_2016.cfg',
    'Emb_16' : 'h-tautau/Production/data/triggers_2016_emb.cfg',
    'MC_17'   : 'h-tautau/Production/data/triggers_2017.cfg',
    'Run2017'   : 'h-tautau/Production/data/triggers_2017.cfg',
    'Emb_17' : 'h-tautau/Production/data/triggers_2017_emb.cfg',
    'MC_18'   : 'h-tautau/Production/data/triggers_2018.cfg',
    'Run2018ABC'   : 'h-tautau/Production/data/triggers_2018.cfg',
    'Run2018D'   : 'h-tautau/Production/data/triggers_2018.cfg',
    'Emb_18ABC' : 'h-tautau/Production/data/triggers_2018_emb.cfg',
    'Emb_18D' : 'h-tautau/Production/data/triggers_2018_emb.cfg'
}

def IsData(sampleType):
    isData = sampleType in dataSampleTypes
    if not isData and not sampleType in mcSampleTypes:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return isData

def GetPeriod(sampleType):
    if sampleType not in periodDict:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return periodDict[sampleType]

def GetGlobalTag(sampleType):
    if sampleType not in globalTagMap:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return globalTagMap[sampleType]

def GetTriggerCfg(sampleType):
    if sampleType not in hltPaths :
        print "ERROR: no HLT paths found for period '{}' .".format(sampleType)
        sys.exit(1)
    return hltPaths[sampleType]
