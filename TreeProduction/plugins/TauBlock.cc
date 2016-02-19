/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TauPFSpecific.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Utilities/General/interface/FileInPath.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "h-tautau/TreeProduction/interface/Tau.h"
#include "h-tautau/TreeProduction/interface/TriggerTools.h"

#define SIMPLE_VAR(type, name) tauTree.name() = patTau.tauID(#name);

class TauBlock : public edm::EDAnalyzer {
public:
    explicit TauBlock(const edm::ParameterSet& iConfig) :
        _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
        tauTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { tauTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    edm::InputTag _inputTag;
    ntuple::TauTree tauTree;
};

void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    tauTree.RunId() = iEvent.id().run();
    tauTree.LumiBlock() = iEvent.id().luminosityBlock();
    tauTree.EventId() = iEvent.id().event();

    edm::Handle<pat::TauCollection> tausHandle;
    iEvent.getByLabel(_inputTag, tausHandle);
    const pat::TauCollection* taus = tausHandle.product();
    if(!taus)
        throw std::runtime_error("Tau collection not found.");

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();

    for (const pat::Tau& patTau : *taus) {
        // Store Tau variables
        tauTree.eta()    = patTau.eta();
        tauTree.phi()    = patTau.phi();
        tauTree.pt()     = patTau.pt();
        tauTree.mass() = patTau.mass();
        tauTree.charge() = patTau.charge();
        tauTree.decayMode() = patTau.decayMode();

        // Leading particle pT
        tauTree.leadChargedParticlePt() = patTau.leadPFChargedHadrCand().isNonnull()
                                          ? patTau.leadPFChargedHadrCand()->pt() : 0.;
        tauTree.leadNeutralParticleEt() = patTau.leadPFNeutralCand().isNonnull()
                                          ? patTau.leadPFNeutralCand()->et() : 0.;
        tauTree.leadParticleEt()        = patTau.leadPFCand().isNonnull()
                                          ? patTau.leadPFCand()->et() : 0.;

        // Number of charged/neutral candidates and photons in different cones
        tauTree.numChargedHadronsSignalCone() = patTau.signalPFChargedHadrCands().size();
        tauTree.numNeutralHadronsSignalCone() = patTau.signalPFNeutrHadrCands().size();
        tauTree.numPhotonsSignalCone()        = patTau.signalPFGammaCands().size();
        tauTree.numParticlesSignalCone()      = patTau.signalPFCands().size();

        tauTree.numChargedHadronsIsoCone() = patTau.isolationPFChargedHadrCands().size();
        tauTree.numNeutralHadronsIsoCone() = patTau.isolationPFNeutrHadrCands().size();
        tauTree.numPhotonsIsoCone()        = patTau.isolationPFGammaCands().size();
        tauTree.numParticlesIsoCone()      = patTau.isolationPFCands().size();

        tauTree.ptSumPFChargedHadronsIsoCone() = patTau.isolationPFChargedHadrCandsPtSum();
        tauTree.etSumPhotonsIsoCone()          = patTau.isolationPFGammaCandsEtSum();

        // tau id. discriminators
        // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
        // for discriminator names see PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py
        TAU_DISCRIMINATOR_DATA()

        tauTree.leadPFCand_mva_e_pi() = patTau.leadPFCand().isNonnull() ? patTau.leadPFCand()->mva_e_pi() : 1.;

        // NEW quantities
        tauTree.emFraction()              = patTau.emFraction();
        tauTree.maximumHCALPFClusterEt()  = patTau.maximumHCALPFClusterEt();
        tauTree.ecalStripSumEOverPLead()  = patTau.ecalStripSumEOverPLead();
        tauTree.bremsRecoveryEOverPLead() = patTau.bremsRecoveryEOverPLead();
        tauTree.hcalTotOverPLead()        = patTau.hcalTotOverPLead();
        tauTree.hcalMaxOverPLead()        = patTau.hcalMaxOverPLead();
        tauTree.hcal3x3OverPLead()        = patTau.hcal3x3OverPLead();

        tauTree.etaetaMoment() = patTau.etaetaMoment();
        tauTree.phiphiMoment() = patTau.phiphiMoment();
        tauTree.etaphiMoment() = patTau.etaphiMoment();

        // Vertex information
        const reco::Candidate::Point& vertex = patTau.vertex();
        tauTree.vx() = vertex.x();
        tauTree.vy() = vertex.y();
        tauTree.vz() = vertex.z();

        // Charged Hadrons Candidates
        for(size_t n = 0; n < patTau.signalPFChargedHadrCands().size(); ++n) {
            tauTree.signalChHadCand_Pt().push_back(patTau.signalPFChargedHadrCands()[n]->pt());
            tauTree.signalChHadCand_Eta().push_back(patTau.signalPFChargedHadrCands()[n]->eta());
            tauTree.signalChHadCand_Phi().push_back(patTau.signalPFChargedHadrCands()[n]->phi());
        }

        for(size_t n = 0; n < patTau.isolationPFChargedHadrCands().size(); ++n) {
            tauTree.isoChHadCand_Pt().push_back(patTau.isolationPFChargedHadrCands()[n]->pt());
            tauTree.isoChHadCand_Eta().push_back(patTau.isolationPFChargedHadrCands()[n]->eta());
            tauTree.isoChHadCand_Phi().push_back(patTau.isolationPFChargedHadrCands()[n]->phi());
        }

        // Neutral Hadrons Candidates
        for(size_t n = 0; n < patTau.signalPFNeutrHadrCands().size(); ++n) {
            tauTree.signalNeutrHadCand_Pt().push_back(patTau.signalPFNeutrHadrCands()[n]->pt());
            tauTree.signalNeutrHadCand_Eta().push_back(patTau.signalPFNeutrHadrCands()[n]->eta());
            tauTree.signalNeutrHadCand_Phi().push_back(patTau.signalPFNeutrHadrCands()[n]->phi());
        }

        for(size_t n = 0; n < patTau.isolationPFNeutrHadrCands().size(); ++n) {
            tauTree.isoNeutrHadCand_Pt().push_back(patTau.isolationPFNeutrHadrCands()[n]->pt());
            tauTree.isoNeutrHadCand_Eta().push_back(patTau.isolationPFNeutrHadrCands()[n]->eta());
            tauTree.isoNeutrHadCand_Phi().push_back(patTau.isolationPFNeutrHadrCands()[n]->phi());
        }

        // Gamma Candidates
        for(size_t n = 0; n < patTau.signalPFGammaCands().size(); ++n) {
            tauTree.signalGammaCand_Pt().push_back(patTau.signalPFGammaCands()[n]->pt());
            tauTree.signalGammaCand_Eta().push_back(patTau.signalPFGammaCands()[n]->eta());
            tauTree.signalGammaCand_Phi().push_back(patTau.signalPFGammaCands()[n]->phi());
        }

        for(size_t n = 0; n < patTau.isolationPFGammaCands().size(); ++n) {
            tauTree.isoGammaCand_Pt().push_back(patTau.isolationPFGammaCands()[n]->pt());
            tauTree.isoGammaCand_Eta().push_back(patTau.isolationPFGammaCands()[n]->eta());
            tauTree.isoGammaCand_Phi().push_back(patTau.isolationPFGammaCands()[n]->phi());
        }

        tauTree.matchedTriggerPaths() = CollectMatchedTriggerPaths(patTau);
        tauTree.Fill();
    }
}

#undef SIMPLE_VAR
#undef TAU_DISCRIMINATOR_DATA

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
