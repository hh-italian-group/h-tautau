# Configurations dependent on the sample type.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import sys
from sets import Set
import FWCore.ParameterSet.Config as cms

mcSampleTypes = Set([ 'Summer16MC', 'Fall17MC' ])
dataSampleTypes = Set([ 'Run2016' , 'Run2017' ])

periodDict = { 'Summer16MC' : 'Run2016',
               'Run2016' : 'Run2016',
               'Fall17MC' : 'Run2017',
               'Run2017' : 'Run2017'
               }

hltPaths = {
    'Run2016'   : 'h-tautau/Production/data/triggers_2016.cfg',
    'Run2017'   : 'h-tautau/Production/data/triggers_2017.cfg'
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

def GetTriggerCfg(period):
    if period not in hltPaths :
        print "ERROR: no HLT paths found for period '{}' .".format(period)
        sys.exit(1)
    return hltPaths[period]



