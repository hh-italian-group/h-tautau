# Produce EventTuple for all channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import importlib
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *

tauIdDiscrMVA_trainings_run2_2017 = {
    'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT2017",
}
tauIdDiscrMVA_WPs_run2_2017 = {
    'tauIdMVAIsoDBoldDMwLT2017' : {
    'Eff95' : "DBoldDMwLTEff95",
    'Eff90' : "DBoldDMwLTEff90",
    'Eff80' : "DBoldDMwLTEff80",
    'Eff70' : "DBoldDMwLTEff70",
    'Eff60' : "DBoldDMwLTEff60",
    'Eff50' : "DBoldDMwLTEff50",
    'Eff40' : "DBoldDMwLTEff40"
}
}
tauIdDiscrMVA_2017_version = "v1"

def loadMVA_WPs_run2_2017(process):
    
    for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017.items():
        
        process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
           cms.PSet(
              record = cms.string('GBRWrapperRcd'),
              tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version)),
              label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version))
           )
        )
            
        for WP in tauIdDiscrMVA_WPs_run2_2017[training].keys():
           process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
              cms.PSet(
                 record = cms.string('PhysicsTGraphPayloadRcd'),
                 tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP)),
                 label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP))
              )
           )
                    
           process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
              cms.PSet(
                 record = cms.string('PhysicsTFormulaPayloadRcd'),
                 tag = cms.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version)),
                 label = cms.untracked.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version))
              )
           )

def tauId_calculation_2017(process):
    
    process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')

    loadMVA_WPs_run2_2017(process)

    process.rerunDiscriminationByIsolationMVArun2v1raw = patDiscriminationByIsolationMVArun2v1raw.clone(
       PATTauProducer = cms.InputTag('slimmedTaus'),
       Prediscriminants = noPrediscriminants,
       loadMVAfromDB = cms.bool(True),
       mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1"),  # name of the training you want to use: training with 2017 MC_v1 for oldDM
       mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
       requireDecayMode = cms.bool(True),
       verbosity = cms.int32(0)
    )

    process.rerunDiscriminationByIsolationMVArun2v1VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
       PATTauProducer = cms.InputTag('slimmedTaus'),
       Prediscriminants = noPrediscriminants,
       toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
       key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw:category'),
       loadMVAfromDB = cms.bool(True),
       mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization"), # normalization fo the training you want to use
       mapping = cms.VPSet(
          cms.PSet(
             category = cms.uint32(0),
             cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90"), # this is the name of the working point you want to use
             variable = cms.string("pt"),
          )
       )
     )

# here we produce all the other working points for the training
    process.rerunDiscriminationByIsolationMVArun2v1VVLoose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95")
    process.rerunDiscriminationByIsolationMVArun2v1Loose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80")
    process.rerunDiscriminationByIsolationMVArun2v1Medium = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70")
    process.rerunDiscriminationByIsolationMVArun2v1Tight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60")
    process.rerunDiscriminationByIsolationMVArun2v1VTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50")
    process.rerunDiscriminationByIsolationMVArun2v1VVTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
    process.rerunDiscriminationByIsolationMVArun2v1VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40")

# this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.
    process.rerunMvaIsolation2SeqRun2 = cms.Sequence(
       process.rerunDiscriminationByIsolationMVArun2v1raw
       *process.rerunDiscriminationByIsolationMVArun2v1VLoose
       *process.rerunDiscriminationByIsolationMVArun2v1VVLoose
       *process.rerunDiscriminationByIsolationMVArun2v1Loose
       *process.rerunDiscriminationByIsolationMVArun2v1Medium
       *process.rerunDiscriminationByIsolationMVArun2v1Tight
       *process.rerunDiscriminationByIsolationMVArun2v1VTight
       *process.rerunDiscriminationByIsolationMVArun2v1VVTight
    )

# embed new id's into new tau collection
    embedID = cms.EDProducer("PATTauIDEmbedder",
       src = cms.InputTag('slimmedTaus'),
       tauIDSources = cms.PSet(
          byIsolationMVArun2v1DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
          byVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VLoose'),
          byVVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVLoose'),
          byLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Loose'),
          byMediumIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Medium'),
          byTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Tight'),
          byVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VTight'),
          byVVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVTight'),
       ),
    )

    setattr(process, "NewTauIDsEmbedded", embedID)

# inclusion in the process
#process.p += process.rerunMvaIsolation2SeqRun2
#process.p += getattr(process, "NewTauIDsEmbedded")
#  then you can continue with your ntuple creation process for example

