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

DeltaSafetyPt = { 'e': 3, 'mu' : 2, 'tau': 5 }

hltPaths_eTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_Ele23_WPLoose_Gsf_v"),
              legs = cms.untracked.vstring('type=e pt=23 filters="hltEle23WPLooseGsfTrackIsoFilter hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"') ),
    cms.PSet( pattern = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
              legs = cms.untracked.vstring('type=e pt=22 filters="hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"' ,
              'type=tau pt=20 filters="hltPFTau20TrackLooseIso hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"') ),
)

hltPaths_eTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_Ele32_WPTight_Gsf_v"),
              legs = cms.untracked.vstring('type=e pt=32 filters="hltEle32WPTightGsfTrackIsoFilter"') ),
    cms.PSet( pattern = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
              legs = cms.untracked.vstring('type=e pt=24 filters="hltEle24erWPTightGsfTrackIsoFilterForTau hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"' ,
                                           'type=tau pt=30 filters="hltSelectedPFTau30LooseChargedIsolationL1HLTMatched hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"') ),
)

hltPaths_eTau = {
    'Run2016'   : hltPaths_eTau_Run2016,
    'Run2017'   : hltPaths_eTau_Run2017
}

hltPaths_muTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu18_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu22_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu18_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu20_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu22_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu17LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu19LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20") ),
)

hltPaths_muTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                               "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27MediumChargedIsolationAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27MediumChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27TightChargedIsolationAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27TightChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_eta2p1_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseChargedIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu24LooseChargedIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseChargedIsoTightOOSCPhotonsAgainstMuon",
                                               "hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_eta2p1_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24MediumChargedIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackMediumChargedIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu24MediumChargedIsoPFTau20") ),
    cms.PSet( pattern = cms.string(" HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackMediumChargedIsoTightOOSCPhotonsAgainstMuon",
                                               "hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_eta2p1_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24TightChargedIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackTightChargedIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu24TightChargedIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                              "hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackTightChargedIsoTightOOSCPhotonsAgainstMuon",
                                               "hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso",
                                               "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
              nLegs = cms.untracked.uint32(1),
              filters2 = cms.untracked.vstring("hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIso1Prong",
                                               "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong") ),
)


hltPaths_muTau = {
    'Run2016'   : hltPaths_muTau_Run2016,
    'Run2017'   : hltPaths_muTau_Run2017
}

hltPaths_tauTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2)),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2)),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v"),
              nLegs = cms.untracked.uint32(2)),
)

hltPaths_tauTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso",
                                               "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched") ),
    cms.PSet( pattern = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong",
                                               "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
              filters2 = cms.untracked.vstring("hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIso1Prong",
                                               "hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong") ),
    cms.PSet( pattern = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau20TrackPt1LooseChargedIsolationReg") ),
    cms.PSet( pattern = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltDoublePFTau20TrackPt1MediumChargedIsolationReg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau20TrackPt1LooseChargedIsolationReg") ),
    cms.PSet( pattern = cms.string("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong"),
              filters2 = cms.untracked.vstring("hltDoublePFTau20TrackPt1LooseChargedIsolationReg") ),

)

hltPaths_tauTau = {
    'Run2016'   : hltPaths_tauTau_Run2016,
    'Run2017'   : hltPaths_tauTau_Run2017
}

hltPaths_muMu_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") )
)

hltPaths_muMu_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07") ),
)

hltPaths_muMu = {
    'Run2016'   : hltPaths_muMu_Run2016,
    'Run2017'   : hltPaths_muMu_Run2017
}

hltPaths = { 'eTau' : hltPaths_eTau, 'muTau' : hltPaths_muTau, 'tauTau' : hltPaths_tauTau, 'muMu' : hltPaths_muMu }
for channelName in hltPaths:
    hltPaths[channelName]['Summer16MC'] = hltPaths[channelName]['Run2016']
    hltPaths[channelName]['Fall17MC'] = hltPaths[channelName]['Run2017']

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

def GetPeriod(sampleType):
    if sampleType not in periodDict:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return periodDict[sampleType]
