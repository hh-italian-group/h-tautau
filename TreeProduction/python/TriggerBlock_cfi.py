## Configuration file that defines the producer of ROOT-tuple for trigger information.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

triggerBlock = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.int32(0),
  l1InputTag  = cms.InputTag('gtDigis'),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  hltPathsOfInterest = cms.vstring ("HLT_DoubleMu",
                                    "HLT_Mu",
                                    "HLT_IsoMu", 
                                    "HLT_TripleMu",
                                    "IsoPFTau",
                                    "TrkIsoT", 
                                    "HLT_Ele")
)
