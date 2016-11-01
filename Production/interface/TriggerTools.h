/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {

enum class CMSSW_Process { SIM, HLT, RECO, PAT };
ENUM_NAMES(CMSSW_Process) = {
    { CMSSW_Process::SIM, "SIM" },
    { CMSSW_Process::HLT, "HLT" },
    { CMSSW_Process::RECO, "RECO" },
    { CMSSW_Process::PAT, "PAT" }
};

class TriggerTools {
public:
    using L1ParticlePtrSet = std::set<const l1extra::L1JetParticle*>;
    template<typename T> using EDGetTokenT = edm::EDGetTokenT<T>;
    template<typename T> using Handle = edm::Handle<T>;

    TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                 EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                 EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                 EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token);

    TriggerTools(const edm::ParameterSet& iConfig);

    void Initialize(const edm::Event& iEvent);

    bool HaveTriggerFired(const std::vector<std::string> &hltPaths);
    const pat::TriggerObjectStandAlone* FindMatchingTriggerObject(const std::string& pathOfInterest,
                                                                  const LorentzVector& candidateMomentum,
                                                                  double deltaR_Limit);

    template<typename HiggsCandidate>
    std::vector<HiggsCandidate> ApplyTriggerMatch(const std::vector<HiggsCandidate>& higgses,
                                                  const std::vector<std::string>& hltPaths, bool isCrossTrigger)
    {
        std::vector<HiggsCandidate> triggeredHiggses;
        for (const auto& higgs : higgses) {
            for (const std::string& hltPath : hltPaths) {
                if(HaveTriggerMatched(hltPath, higgs, cuts::H_tautau_2016::DeltaR_triggerMatch, isCrossTrigger)) {
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

    bool TryGetTriggerResult(CMSSW_Process process, const std::string& name, bool& result) const;
    bool GetTriggerResult(CMSSW_Process process, const std::string& name) const;
    bool GetAnyTriggerResult(const std::string& name) const;

private:
    std::map<CMSSW_Process, EDGetTokenT<edm::TriggerResults>> triggerResults_tokens;
    EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_token;
    EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_token;
    EDGetTokenT<std::vector<l1extra::L1JetParticle>> l1JetParticles_token;

    const edm::Event* iEvent;
    std::map<CMSSW_Process, Handle<edm::TriggerResults>> triggerResultsMap;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;
};

} // namespace analysis
