# Configurations dependent on the sample type.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import sys
from sets import Set

mcSampleTypes = Set([ 'Spring15MC', 'Fall15MC', 'Spring16MC' ])
dataSampleTypes = Set([ 'Run2015B', 'Run2015C', 'Run2015D' , 'Run2016B'])

hltPaths_eTau = {
    'Spring15MC' : [ "HLT_Ele22_eta2p1_WP75_Gsf" ],
    'Fall15MC'   : [ "HLT_Ele23_WPLoose_Gsf_v3" ],
    'Spring16MC' : [ "HLT_Ele25_eta2p1_WPTight_Gsf_v2" ],
    'Run2015C'   : [ "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1", "HLT_Ele32_eta2p1_WPTight_Gsf_v1" ],
    'Run2015D'   : [ "HLT_Ele23_WPLoose_Gsf_v" ],
    'Run2016B'   : [ "HLT_Ele25_eta2p1_WPTight_Gsf_v" ]
}

hltPaths_muTau = {
    'Spring15MC' : [ "HLT_IsoMu17_eta2p1_v1" ],
    'Fall15MC'   : [ "HLT_IsoMu18_v2" ],
    'Spring16MC' : [ "HLT_IsoMu22_v3" ],
    'Run2015B'   : [ "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2", "HLT_IsoMu24_eta2p1_v2" ],
    'Run2015C'   : [ "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2", "HLT_IsoMu24_eta2p1_v2" ],
    'Run2015D'   : [ "HLT_IsoMu18_v" ],
    'Run2016B'   : [ "HLT_IsoMu22_v" ]
}

hltPaths_tauTau = {
    'Spring16MC' : [ "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v3" ],
    'Run2015C' : [ "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2" ],
    'Run2015D' : [ "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v" ],
    'Run2016B' : [ "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v" ]
}

hltPaths = { 'eTau' : hltPaths_eTau, 'muTau' : hltPaths_muTau, 'tauTau' : hltPaths_tauTau }

def IsData(sampleType):
    isData = sampleType in dataSampleTypes
    if not isData and not sampleType in mcSampleTypes:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return isData

def GetHltPaths(channelName, sampleType):
    if not channelName in hltPaths or not sampleType in hltPaths[channelName]:
        print "ERROR: no HLT paths found for sample type '{}' for channel '{}'.".format(sampleType, channel)
        sys.exit(1)
    return hltPaths[channelName][sampleType]
