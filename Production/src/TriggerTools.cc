/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TriggerTools.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


namespace analysis {

TriggerTools::TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResults_token,
                           EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                           EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                           EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token) :
    triggerResults_token(_triggerResults_token), triggerPrescales_token(_triggerPrescales_token),
    triggerObjects_token(_triggerObjects_token), l1JetParticles_token(_l1JetParticles_token)
{
}


//TriggerTools::TriggerTools(const edm::ParameterSet& iConfig) :
//    triggerResults_token(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
//    triggerPrescales_token(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
//    triggerObjects_token(consumes<pat::TriggerObjectStandAloneCollection>(
//                            iConfig.getParameter<edm::InputTag>("objects"))),
//    l1JetParticles_token(mayConsume<std::vector<l1extra::L1JetParticle>>(
//                            iConfig.getParameter<edm::InputTag>("l1JetParticleProduct")))
//{

//}

void TriggerTools::Initialize(const edm::Event &_iEvent)
{
    iEvent = &_iEvent;
    iEvent->getByToken(triggerResults_token, triggerResults);
    iEvent->getByToken(triggerPrescales_token, triggerPrescales);
    iEvent->getByToken(triggerObjects_token, triggerObjects);
    iEvent->getByToken(l1JetParticles_token, l1JetParticles);
}

bool TriggerTools::HaveTriggerFired(const std::set<std::string>& hltPaths)
{
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResults);

    for (unsigned i = 0; i < triggerResults->size(); ++i) {
        if(triggerPrescales->getPrescaleForIndex(i) != 1 && triggerResults->accept(i)) continue;
        const std::string& objectMatchedPath = triggerNames.triggerName(i);
        for (const std::string& triggerPath : hltPaths ){
            const size_t found = objectMatchedPath.find(triggerPath);
            if(found != std::string::npos)
                return true;
        }
    }
    return false;
}

const pat::TriggerObjectStandAlone* TriggerTools::FindMatchingTriggerObject(const std::string& pathOfInterest,
                                                                            const LorentzVector& candidateMomentum,
                                                                            double deltaR_Limit)
{
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResults);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        if(ROOT::Math::VectorUtil::DeltaR2(triggerObject.polarP4(), candidateMomentum) >= deltaR2) continue;
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        const auto& matchedPaths = unpackedTriggerObject.pathNames(true, true);
        for(const auto& matchedPath : matchedPaths) {
            if(matchedPath.find(pathOfInterest) != std::string::npos)
                return &triggerObject;
        }
    }
    return nullptr;
}

TriggerTools::L1ParticlePtrSet TriggerTools::L1TauMatch(const LorentzVector& tauMomentum)
{
    static constexpr double MinPt = 28;
    static constexpr double DeltaR = 0.5;
    static constexpr double DeltaR2 = std::pow(DeltaR, 2);

    L1ParticlePtrSet result;
    for(const auto& l1Tau : *l1JetParticles){
        if(l1Tau.pt() > MinPt && ROOT::Math::VectorUtil::DeltaR2(l1Tau.polarP4(), tauMomentum) < DeltaR2)
            result.insert(&l1Tau);
    }
    return result;
}

} // namespace analysis
