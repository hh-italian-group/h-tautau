/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"

namespace analysis {

class TriggerTools {
public:
    using L1ParticlePtrSet = std::set<const l1extra::L1JetParticle*>;
    template<typename T> using EDGetTokenT = edm::EDGetTokenT<T>;

    TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResults_token,
                 EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                 EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                 EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token);

    TriggerTools(const edm::ParameterSet& iConfig);

    void Initialize(const edm::Event& iEvent);

    bool HaveTriggerFired(const std::set<std::string>& hltPaths);
    const pat::TriggerObjectStandAlone* FindMatchingTriggerObject(const std::string& pathOfInterest,
                                                                  const LorentzVector& candidateMomentum,
                                                                  double deltaR_Limit);

    template<typename HiggsCandidate>
    std::vector<HiggsCandidate> ApplyTriggerMatch(const std::vector<HiggsCandidate>& higgses,
                                                  const std::set<std::string>& hltPaths, bool isCrossTrigger)
    {
        std::vector<HiggsCandidate> triggeredHiggses;
        for (const auto& higgs : higgses) {
            for (const std::string& hltPath : hltPaths) {
                if(HaveTriggerMatched(hltPath, higgs, cuts::Htautau_2015::DeltaR_triggerMatch, isCrossTrigger)) {
                    triggeredHiggses.push_back(higgs);
                    break;
                }
            }
        }
        return triggeredHiggses;
    }

    template<typename HiggsCandidate>
    bool HaveTriggerMatched(const std::string& pathOfInterest, const HiggsCandidate& candidate,
                            double deltaR_Limit, bool isCrossTrigger)
    {
        if(!FindMatchingTriggerObject(pathOfInterest, candidate.GetFirstDaughter().GetMomentum(), deltaR_Limit))
            return false;
        return !isCrossTrigger
                || FindMatchingTriggerObject(pathOfInterest, candidate.GetSecondDaughter().GetMomentum(), deltaR_Limit);
    }

    L1ParticlePtrSet L1TauMatch(const LorentzVector& tauMomentum);

    template<typename HiggsCandidate>
    std::vector<HiggsCandidate> ApplyL1TriggerTauMatch(const std::vector<HiggsCandidate>& higgses)
    {
        std::vector<HiggsCandidate> L1Higgses;
        for (const auto& higgs : higgses) {
            const auto leg1_match = L1TauMatch(higgs.GetFirstDaughter().GetMomentum());
            if(!leg1_match.size()) continue;
            const auto leg2_match = L1TauMatch(higgs.GetSecondDaughter().GetMomentum());
            if(!leg2_match.size()) continue;

            if(leg1_match.size() > 1 || leg2_match.size() > 1 || *leg1_match.begin() != *leg2_match.begin())
                L1Higgses.push_back(higgs);
        }
        return L1Higgses;
    }

private:
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_token;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_token;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_token;
    edm::EDGetTokenT<std::vector<l1extra::L1JetParticle>> l1JetParticles_token;

    const edm::Event* iEvent;
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;
};

} // namespace analysis
