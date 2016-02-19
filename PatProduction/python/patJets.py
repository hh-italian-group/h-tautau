##  Configuration file that defines parameters related to PAT Jet objects.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.jetTools import *
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

def applyJetParameters(process, isMC, isEmbedded):
    # Jet corrections
    jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
    if not isMC:
        jec.extend([ 'L2L3Residual' ])
    switchJetCollection(process, cms.InputTag('ak5PFJets'),
         doJTA        = True,
         doBTagging   = True,
         jetCorrLabel = ('AK5PF', cms.vstring(jec)),
         doType1MET   = True,
         genJetCollection=cms.InputTag("ak5GenJets"),
         doJetID      = True,
         jetIdLabel   = 'ak5'
    )

    if isEmbedded:
        process.jetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")

    process.load("RecoJets.JetProducers.PileupJetID_cfi")

    process.puJetId = cms.EDProducer("PileupJetIdProducer",
        residualsTxt = cms.FileInPath('HHbbTauTau/PatProduction/data/dummy.txt'),
        runMvas = cms.bool(False),
        inputIsCorrected = cms.bool(False),
        vertexes = cms.InputTag("offlinePrimaryVertices"),
        produceJetIds = cms.bool(True),
        jec = cms.string('AK5PF'),
        residualsFromTxt = cms.bool(False),
        applyJec = cms.bool(True),
        jetids = cms.InputTag(""),
        rho = cms.InputTag("kt6PFJets","rho"),
        jets = cms.InputTag("ak5PFJets"),
        algos = cms.VPSet(cms.PSet(
            cutBased = cms.bool(True),
            JetIdParams = cms.PSet(
                Pt010_RMSLoose = cms.vdouble(0.06, 0.05, 0.05, 0.07),
                Pt3050_RMSMedium = cms.vdouble(0.06, 0.03, 0.03, 0.04),
                Pt3050_BetaStarMedium = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt010_BetaStarMedium = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt1020_BetaStarTight = cms.vdouble(0.15, 0.15, 999.0, 999.0),
                Pt010_RMSMedium = cms.vdouble(0.06, 0.03, 0.03, 0.04),
                Pt1020_RMSLoose = cms.vdouble(0.06, 0.05, 0.05, 0.07),
                Pt1020_BetaStarLoose = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt2030_RMSTight = cms.vdouble(0.05, 0.07, 0.03, 0.045),
                Pt3050_RMSTight = cms.vdouble(0.05, 0.06, 0.03, 0.04),
                Pt1020_RMSTight = cms.vdouble(0.06, 0.07, 0.04, 0.05),
                Pt3050_BetaStarTight = cms.vdouble(0.15, 0.15, 999.0, 999.0),
                Pt3050_RMSLoose = cms.vdouble(0.06, 0.05, 0.05, 0.055),
                Pt2030_RMSLoose = cms.vdouble(0.06, 0.05, 0.05, 0.055),
                Pt010_BetaStarTight = cms.vdouble(0.15, 0.15, 999.0, 999.0),
                Pt2030_BetaStarMedium = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt1020_RMSMedium = cms.vdouble(0.06, 0.03, 0.03, 0.04),
                Pt2030_BetaStarLoose = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt2030_BetaStarTight = cms.vdouble(0.15, 0.15, 999.0, 999.0),
                Pt2030_RMSMedium = cms.vdouble(0.06, 0.03, 0.03, 0.04),
                Pt010_BetaStarLoose = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt3050_BetaStarLoose = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt1020_BetaStarMedium = cms.vdouble(0.2, 0.3, 999.0, 999.0),
                Pt010_RMSTight = cms.vdouble(0.06, 0.07, 0.04, 0.05)
            ),
            impactParTkThreshold = cms.double(1.0),
            label = cms.string('cutbased')
        ))
    )

    process.pileupJetIdProducer = cms.EDProducer("PileupJetIdProducer",
        produceJetIds = cms.bool(False),
        runMvas = cms.bool(True),
        inputIsCorrected = cms.bool(False),
        vertexes = cms.InputTag("offlinePrimaryVertices"),
        residualsTxt = cms.FileInPath('HHbbTauTau/PatProduction/data/dummy.txt'),
        jec = cms.string('AK5PF'),
        residualsFromTxt = cms.bool(False),
        applyJec = cms.bool(True),
        jetids = cms.InputTag("puJetId"),
        rho = cms.InputTag("kt6PFJets","rho"),
        jets = cms.InputTag("ak5PFJets"),
        algos = cms.VPSet(cms.PSet(
            tmvaVariables = cms.vstring('nvtx',
                'dZ',
                'beta',
                'betaStar',
                'nCharged',
                'nNeutrals',
                'dR2Mean',
                'ptD',
                'frac01',
                'frac02',
                'frac03',
                'frac04',
                'frac05'),
            tmvaMethod = cms.string('JetIDMVAHighPt'),
            cutBased = cms.bool(False),
            tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),
            tmvaSpectators = cms.vstring('jetPt',
                'jetEta',
                'jetPhi'),
            label = cms.string('full'),
            version = cms.int32(-1),
            JetIdParams = cms.PSet(
                Pt2030_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42),
                Pt2030_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
                Pt3050_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
                Pt1020_MET = cms.vdouble(0.3, -0.2, -0.4, -0.4),
                Pt2030_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
                Pt010_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
                Pt1020_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
                Pt3050_MET = cms.vdouble(0.0, 0.0, -0.1, -0.2),
                Pt010_MET = cms.vdouble(0.0, -0.6, -0.4, -0.4),
                Pt1020_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
                Pt010_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
                Pt1020_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
                Pt2030_MET = cms.vdouble(0.0, 0.0, 0.0, 0.0),
                Pt010_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
                Pt3050_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
                Pt3050_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42)
            ),
            impactParTkThreshold = cms.double(1.0)
        ),
            cms.PSet(
                tmvaVariables = cms.vstring('nvtx',
                    'jetPt',
                    'jetEta',
                    'jetPhi',
                    'dZ',
                    'beta',
                    'betaStar',
                    'nCharged',
                    'nNeutrals',
                    'dR2Mean',
                    'ptD',
                    'frac01',
                    'frac02',
                    'frac03',
                    'frac04',
                    'frac05'),
                tmvaMethod = cms.string('JetID'),
                cutBased = cms.bool(False),
                tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml'),
                tmvaSpectators = cms.vstring(),
                label = cms.string('met'),
                version = cms.int32(-1),
                JetIdParams = JetIdParams,
                impactParTkThreshold = cms.double(0.)
            )
        )
    )

    # Embed into PAT jets as userdata
    process.patJets.userData.userFloats.src = cms.VInputTag(
        cms.InputTag('pileupJetIdProducer', 'fullDiscriminant'),
        cms.InputTag('pileupJetIdProducer', 'metDiscriminant')
    )
    process.patJets.userData.userInts.src = cms.VInputTag(
        cms.InputTag('pileupJetIdProducer', 'fullId'),
        cms.InputTag('pileupJetIdProducer', 'metId')
    )

    process.PatJetsWithEmbeddedVariables = cms.EDProducer('JetUserEmbedder',
        jetTag = cms.InputTag("patJetsTriggerMatch"),
        corrector = cms.string("ak5PFL1Fastjet")
    )

    return
