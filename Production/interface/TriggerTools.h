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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "h-tautau/Analysis/include/TriggerResults.h"

namespace analysis {

enum class CMSSW_Process { SIM, HLT, RECO, PAT };
ENUM_NAMES(CMSSW_Process) = {
    { CMSSW_Process::SIM, "SIM" },
    { CMSSW_Process::HLT, "HLT" },
    { CMSSW_Process::RECO, "RECO" },
    { CMSSW_Process::PAT, "PAT" }
};

namespace detail {

template<typename PatObject>
LegType GetTriggerObjectTypes(const PatObject&);

template<>
inline LegType GetTriggerObjectTypes<pat::Electron>(const pat::Electron&)
{
    return LegType::e;
}

template<>
inline LegType GetTriggerObjectTypes<pat::Muon>(const pat::Muon&)
{
    return LegType::mu;
}

template<>
inline LegType GetTriggerObjectTypes<pat::Tau>(const pat::Tau&)
{
    return LegType::tau;
}

} // namespace detail

class TriggerTools {
public:
    using L1ParticlePtrSet = std::set<const l1extra::L1JetParticle*>;
    template<typename T> using EDGetTokenT = edm::EDGetTokenT<T>;
    template<typename T> using Handle = edm::Handle<T>;
    using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;

    TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                 EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                 EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                 EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token,
                 std::string _triggerCfg, analysis::Channel _channel);

    TriggerTools(const edm::ParameterSet& iConfig, const std::string& _triggerCfg, const analysis::Channel& _channel);

    static TriggerDescriptorCollection CreateTriggerDescriptors(const std::map<std::string, std::vector<std::string>>& pattern_legs_map,
                                                                const Channel& channel);

    void Initialize(const edm::Event& iEvent);

    void SetTriggerAcceptBits(analysis::TriggerResults& results);

    std::map<size_t,TriggerTools::TriggerObjectSet> FindMatchingTriggerObjects(
            const analysis::TriggerDescriptorCollection::TriggerDescriptor& pattern_struct,
            const LorentzVector& candidateMomentum, const LegType& candidate_type, double deltaR_Limit);


    template<typename HiggsCandidate>
    void SetTriggerMatchBits(analysis::TriggerResults& results, const HiggsCandidate& candidate, double deltaR_Limit)
    { 
        for (size_t n = 0; n < triggerDescriptors.Size(); ++n){
            analysis::TriggerDescriptorCollection::TriggerDescriptor pattern_struct = triggerDescriptors.GetTriggerDescriptor(n);
            std::map<size_t, TriggerObjectSet> matches_first, matches_second;
            matches_first = FindMatchingTriggerObjects(pattern_struct,candidate.GetFirstDaughter().GetMomentum(),
                                                        detail::GetTriggerObjectTypes(*candidate.GetFirstDaughter()),
                                                        deltaR_Limit);
            matches_second = FindMatchingTriggerObjects(pattern_struct,candidate.GetSecondDaughter().GetMomentum(),
                                                         detail::GetTriggerObjectTypes(*candidate.GetSecondDaughter()),
                                                         deltaR_Limit);
            results.SetMatch(n, TriggerMatchFound(matches_first, matches_second, pattern_struct.legs_info.size()));
        }
    }

    //to be fixed
    bool TriggerMatchFound(const std::map<size_t, TriggerObjectSet>& matches_first,
                           const std::map<size_t, TriggerObjectSet>& matches_second,
                           const size_t n_legs_total)
    {
        std::vector<const pat::TriggerObjectStandAlone*> comb_match;
        for(const auto& first : matches_first){
            for(const auto& second : matches_second){
                std::set_difference(first.second.begin(), first.second.end(),
                                      second.second.begin(), second.second.end(),
                                       std::inserter(comb_match,comb_match.begin()));
            }
        }
        if(comb_match.size() >= n_legs_total)
            return true;
        return false;
    }

    bool TryGetTriggerResult(CMSSW_Process process, const std::string& name, bool& result) const;
    bool GetTriggerResult(CMSSW_Process process, const std::string& name) const;
    bool TryGetAnyTriggerResult(const std::string& name, bool& result) const;
    bool GetAnyTriggerResult(const std::string& name) const;

private:
    std::map<CMSSW_Process, EDGetTokenT<edm::TriggerResults>> triggerResults_tokens;
    EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_token;
    EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_token;
    EDGetTokenT<std::vector<l1extra::L1JetParticle>> l1JetParticles_token;

    const edm::Event* iEvent;
    analysis::TriggerDescriptorCollection triggerDescriptors;
    std::string triggerCfg;
    analysis::Channel channel;
    std::set<analysis::Channel> channels;
    std::map<std::string, std::vector<std::string>> pattern_legs_map;
    std::map<analysis::LegType, double> deltaPt_map;
    std::map<CMSSW_Process, Handle<edm::TriggerResults>> triggerResultsMap;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;

    std::map<std::string,std::map<size_t,TriggerObjectSet>> path_legId_triggerObjPtr_map;
};

} // namespace analysis
