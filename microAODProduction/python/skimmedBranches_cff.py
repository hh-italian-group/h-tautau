import FWCore.ParameterSet.Config as cms

BaseMCBranches = cms.untracked.vstring(['drop *',
                                         'keep GenEventInfoProduct_generator__SIM',
                                         'keep edmTriggerResults_*__HLT',
                                         'keep edmTriggerResults_*__PAT',
                                         'keep PileupSummaryInfos_addPileupInfo__HLT',
                                         'keep PileupSummaryInfos_slimmedAddPileupInfo__PAT',
                                         'keep double_*__RECO',
                                         'keep recoBeamSpot_offlineBeamSpot__RECO',
                                         'keep *_electronMVAValueMapProducer_*_USER', ## Value Map with the electrons ID MVA values
                                         'keep floatedmValueMap_offlineSlimmedPrimaryVertices__PAT',
                                         'keep patPackedTriggerPrescales_patTrigger__PAT',
                                         'keep patElectrons_slimmedElectrons__PAT',
                                         'keep patJets_*__PAT',
                                         'keep patMETs_*__PAT', ## MET
                                         'keep patMuons_slimmedMuons__PAT',
                                         'keep patPackedCandidates_packedPFCandidates__PAT',
                                         'keep patPackedGenParticles_packedGenParticles__PAT',
                                         'keep patTaus_slimmedTaus__PAT',
                                         'keep patTriggerObjectStandAlones_selectedPatTrigger__PAT',
                                         'keep recoGenJets_slimmedGenJets__PAT',
                                         'keep recoGenParticles_prunedGenParticles__PAT',
                                         'keep recoVertexs_offlineSlimmedPrimaryVertices__PAT',
                                         'keep recoVertexCompositePtrCandidates_slimmedSecondaryVertices__PAT',
                                         'keep *_l1extraParticles_*_*',
                                         'keep *_METSignificance_*_USER'])

BaseDATABranches = cms.untracked.vstring(['drop *',
                                          'keep edmTriggerResults_*__HLT',
                                          'keep double_*__RECO',
                                          'keep *_electronMVAValueMapProducer_*_USER', ## Value Map with the electrons ID MVA values
                                          'keep floatedmValueMap_offlineSlimmedPrimaryVertices__RECO',
                                          'keep patPackedTriggerPrescales_patTrigger__RECO',
                                          'keep patElectrons_slimmedElectrons__RECO',
                                          'keep patJets_*__RECO',
                                          'keep patMETs_*__RECO', ## MET
                                          'keep patMuons_slimmedMuons__RECO',
                                          'keep patPackedCandidates_packedPFCandidates__RECO',
                                          'keep patPackedGenParticles_packedGenParticles__RECO',
                                          'keep patTaus_slimmedTaus__RECO',
                                          'keep patTriggerObjectStandAlones_selectedPatTrigger__RECO',
                                          'keep recoGenJets_slimmedGenJets__RECO',
                                          'keep recoGenParticles_prunedGenParticles__RECO',
                                          'keep recoVertexs_offlineSlimmedPrimaryVertices__RECO',
                                          'keep recoVertexCompositePtrCandidates_slimmedSecondaryVertices__RECO',
                                          'keep *_l1extraParticles_*_*',
                                          'keep *_METSignificance_*_USER'])
