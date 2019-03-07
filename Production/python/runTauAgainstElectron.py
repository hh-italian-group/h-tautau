import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PATTauDiscriminationAgainstElectronMVA6_cfi import *
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants

def rerunAgainstElectron(process, tauIDsEmbedded):
    process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
    process.rerunDiscriminationAgainstElectronMVA6 = patTauDiscriminationAgainstElectronMVA6.clone(
       PATTauProducer = cms.InputTag('slimmedTaus'),
       Prediscriminants = noPrediscriminants,
       #Prediscriminants = requireLeadTrack,
       loadMVAfromDB = cms.bool(True),
       returnMVA = cms.bool(True),
       method = cms.string("BDTG"),
       mvaName_NoEleMatch_woGwoGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL"),
       mvaName_NoEleMatch_wGwoGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL"),
       mvaName_woGwGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL"),
       mvaName_wGwGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL"),
       mvaName_NoEleMatch_woGwoGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC"),
       mvaName_NoEleMatch_wGwoGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC"),
       mvaName_woGwGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC"),
       mvaName_wGwGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC"),
       minMVANoEleMatchWOgWOgsfBL = cms.double(0.0),
       minMVANoEleMatchWgWOgsfBL  = cms.double(0.0),
       minMVAWOgWgsfBL            = cms.double(0.0),
       minMVAWgWgsfBL             = cms.double(0.0),
       minMVANoEleMatchWOgWOgsfEC = cms.double(0.0),
       minMVANoEleMatchWgWOgsfEC  = cms.double(0.0),
       minMVAWOgWgsfEC            = cms.double(0.0),
       minMVAWgWgsfEC             = cms.double(0.0),
       srcElectrons = cms.InputTag('slimmedElectrons'),
       usePhiAtEcalEntranceExtrapolation = cms.bool(True)
    )
    #tauIDsEmbedded.tauIDSources.againstElectronMVA6RawNew = cms.InputTag('rerunDiscriminationAgainstElectronMVA6')
    #tauIDsEmbedded.tauIDSources.againstElectronMVA6categoryNew = \
    #    cms.InputTag("rerunDiscriminationAgainstElectronMVA6:category")

    # embedID = cms.EDProducer("PATTauIDEmbedder",
    #     src = cms.InputTag('slimmedTaus'),
    #     tauIDSources = cms.PSet(
    #         againstElectronMVA6RawNew = cms.InputTag('rerunDiscriminationAgainstElectronMVA6'),
    #         againstElectronMVA6categoryNew = cms.InputTag("rerunDiscriminationAgainstElectronMVA6:category")
    #         ),
    #     )
    #self.process.NewTauIDsEmbedded = embedID
    #setattr(process, "NewTauIDsEmbedded", embedID)
