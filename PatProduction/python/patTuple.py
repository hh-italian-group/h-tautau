##  Configuration file to produce PAT-tuples and ROOT-tuples for X->HH->bbTauTau analysis.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

process = cms.Process("muAOD")

## Import skeleton process.
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options.wantSummary = False

## Drop input reco taus to ensure that the new tau collection will be used.
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoPFTaus_*_*_*'
    )
#process.source.eventsToProcess = cms.untracked.VEventRange('190706:11370491-190706:11370491')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('190706:15-190706:15')

from HHbbTauTau.PatProduction.patOptions import *
parseAndApplyOptions(process)


from HHbbTauTau.PatProduction.patMuons import *
applyMuonParameters(process)

from HHbbTauTau.PatProduction.patTaus import *
applyTauParameters(process)

from HHbbTauTau.PatProduction.patElectrons import *
applyElectronParameters(process, options.isMC)

from HHbbTauTau.PatProduction.patJets import *
applyJetParameters(process, options.isMC, options.isEmbedded)

from HHbbTauTau.PatProduction.patMET import *
applyMETParameters(process, options.isMC)

from HHbbTauTau.PatProduction.patVertices import *
applyVertexParameters(process)

from HHbbTauTau.PatProduction.patTrigger import *
applyTriggerParameters(process)

from HHbbTauTau.PatProduction.skimFilter import *
applySkim(process)

if options.runTree:
    from HHbbTauTau.PatProduction.treeProduction import *
    addTreeSequence(process, options.includeSim, options.treeOutput)
else:
    process.mainTreeContentSequence = cms.Sequence()
    process.simTreeContentSequence = cms.Sequence()

## Remove MC matching from the default sequence
if not options.isMC:
    # workaround for ill-written removeMCMatching code...
    process.patJetsTriggerMatch.addGenPartonMatch   = cms.bool(False)
    process.patJetsTriggerMatch.embedGenPartonMatch = cms.bool(False)
    process.patJetsTriggerMatch.genPartonMatch      = cms.InputTag('')
    process.patJetsTriggerMatch.addGenJetMatch      = cms.bool(False)
    process.patJetsTriggerMatch.genJetMatch         = cms.InputTag('')
    process.patJetsTriggerMatch.getJetMCFlavour     = cms.bool(False)
    process.patJetsTriggerMatch.JetPartonMapSource  = cms.InputTag('')

    removeMCMatching(process, ['All'])
    process.patDefaultSequence.remove(process.patJetPartonMatch)
    process.patDefaultSequence.remove(process.patJetFlavourId)
    process.patDefaultSequence.remove(process.patJetPartons)
    process.patDefaultSequence.remove(process.patJetPartonAssociation)
    process.patDefaultSequence.remove(process.patJetFlavourAssociation)
    runOnData(process)


## Define process path.
process.p = cms.Path(
    process.PFTau *
    process.pfParticleSelectionSequence *
    process.muIsoSequence *
    process.electronIsoSequence *
    process.mvaTrigV0 *
    process.mvaNonTrigV0 *
    process.type0PFMEtCorrection *
    process.producePFMETCorrections *
    process.puJetId *
    process.pileupJetIdProducer *
    process.patDefaultSequence *
    process.patMuonsWithEmbeddedVariables *
    process.patElectronsWithEmbeddedVariables *
    process.PatJetsWithEmbeddedVariables *
    process.patVertices *
    process.bbttSkim *
    process.mainTreeContentSequence *
    process.simTreeContentSequence
)

## Output selection.
process.out.outputCommands = [ 'drop *' ]

if options.keepPat:
    process.out.outputCommands.extend([
        'keep patElectrons_patElectronsWithEmbeddedVariables_*_*',
        'keep patJets_PatJetsWithEmbeddedVariables_*_*',
        'keep patMETs_*_*_*',
        'keep patMuons_patMuonsWithEmbeddedVariables_*_*',
        'keep patTaus_patTausTriggerMatch_*_*',
        'keep patVertexs_patVertices_*_*',
        'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep *_TriggerResults_*_HLT',
        'keep L1GlobalTriggerReadoutRecord_*_*_*',
        'keep GenFilterInfo_*_*_*',
        'keep *_patTrigger_*_*',
        'keep patTriggerEvent_*_*_*',
        'keep *_particleFlow_*_*'
    ])
    if options.includeSim:
        process.out.outputCommands.extend([
            'keep recoGenParticles_genParticles_*_*',
            'keep *_genMetTrue_*_*'
        ])
else:
    process.out.dropMetaData = cms.untracked.string('ALL')

#print process.dumpPython()
