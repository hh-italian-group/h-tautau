##  Configuration file that defines parameters related to PAT Trigger objects.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

# Trigger and Trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import *

def applyTriggerParameters(process):
    switchOnTrigger(process)

    ## PAT trigger matching
    ## muons
    process.HLTMatchPatMuons = cms.EDProducer(
        # matching in DeltaR, sorting by best DeltaR
        "PATTriggerMatcherDRLessByR"
        # matcher input collections
        , src     = cms.InputTag( 'patMuons' )
        , matched = cms.InputTag( 'patTrigger' )
        # selections of trigger objects
        , matchedCuts = cms.string( 'type( "TriggerMuon" ) && ' +
                                   '( path( "HLT_IsoMu17_eta2p1_LooseIsoPFTau20*" ) || ' +
                                   'path( "HLT_IsoMu18_eta2p1_LooseIsoPFTau20*" ) )' )
        # selection of matches
        , maxDPtRel   = cms.double( 0.5 ) # no effect here
        , maxDeltaR   = cms.double( 0.5 )
        , maxDeltaEta = cms.double( 0.2 ) # no effect here
        # definition of matcher output
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
      )

    ## electrons
    process.HLTMatchPatElectrons = cms.EDProducer(
        # matching in DeltaR, sorting by best DeltaR
        "PATTriggerMatcherDRLessByR"
        # matcher input collections
        , src     = cms.InputTag( 'patElectrons' )
        , matched = cms.InputTag( 'patTrigger' )
        # selections of trigger objects
        , matchedCuts = cms.string( 'type( "TriggerElectron" ) && ' +
                                 '( path( "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20*" ) || ' +
                                 'path( "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20*" ) )' )
        # selection of matches
        , maxDPtRel   = cms.double( 0.5 ) # no effect here
        , maxDeltaR   = cms.double( 0.5 )
        , maxDeltaEta = cms.double( 0.2 ) # no effect here
        # definition of matcher output
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
    )

    ## taus
    process.HLTMatchPatTaus = cms.EDProducer(
        # matching in DeltaR, sorting by best DeltaR
        "PATTriggerMatcherDRLessByR"
        # matcher input collections
        , src     = cms.InputTag( 'patTaus' )
        , matched = cms.InputTag( 'patTrigger' )
        # selections of trigger objects
        , matchedCuts = cms.string( 'type( "TriggerTau" ) && ' +
                                   '( path( "HLT_IsoMu17_eta2p1_LooseIsoPFTau20*" ) || ' +
                                   'path( "HLT_IsoMu18_eta2p1_LooseIsoPFTau20*" ) || ' +
                                   'path( "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20*" ) || ' +
                                   'path( "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1*" ) )' )
        # selection of matches
        , maxDPtRel   = cms.double( 0.5 ) # no effect here
        , maxDeltaR   = cms.double( 0.5 )
        , maxDeltaEta = cms.double( 0.2 ) # no effect here
        # definition of matcher output
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
      )

    ## jets
    process.HLTMatchPatJets = cms.EDProducer(
        # matching in DeltaR, sorting by best DeltaR
        "PATTriggerMatcherDRLessByR"
        # matcher input collections
        , src     = cms.InputTag( 'patJets' )
        , matched = cms.InputTag( 'patTrigger' )
        # selections of trigger objects
        , matchedCuts = cms.string( 'type( "TriggerJet" ) && ' +
                                   '( path( "HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30*" ) || ' +
                                   'path( "HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30*" ) )')
        # selection of matches
        , maxDPtRel   = cms.double( 0.5 ) # no effect here
        , maxDeltaR   = cms.double( 0.5 )
        , maxDeltaEta = cms.double( 0.2 ) # no effect here
        # definition of matcher output
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
      )

    switchOnTriggerMatchEmbedding(process, triggerMatchers = [
                  'HLTMatchPatMuons',
                  'HLTMatchPatElectrons',
                  'HLTMatchPatTaus',
                  'HLTMatchPatJets'
          ])

    return
