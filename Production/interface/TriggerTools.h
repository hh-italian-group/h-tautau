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

inline const std::map<std::string,std::set<trigger::TriggerObjectType>>& GetTriggerObjectTypes()
{
    static const std::map<analysis::LegType,std::set<trigger::TriggerObjectType>> types
            = {
               { analysis::LegType::e, { trigger::TriggerElectron, trigger::TriggerCluster } },
               { analysis::LegType::mu, { trigger::TriggerMuon } },
               { analysis::LegType::tau, { trigger::TriggerTau } }
              };
    return types;
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
                 std::vector<edm::ParameterSet> _hltPaths);

    TriggerTools(const edm::ParameterSet& iConfig);

    void CreateTriggerDescriptors();

    void Initialize(const edm::Event& iEvent);

    void SetTriggerAcceptBits(const analysis::TriggerDescriptors& descriptors, analysis::TriggerResults& results);

    TriggerObjectSet FindMatchingTriggerObjects(const analysis::TriggerDescriptors& descriptors, size_t path_index,
            const std::set<trigger::TriggerObjectType>& objectTypes, const LorentzVector& candidateMomentum,
            size_t leg_id, double deltaR_Limit);

    template<typename Candidate>
    TriggerObjectSet FindMatchingTriggerObjects(const analysis::TriggerDescriptors& descriptors, size_t path_index,
            const Candidate& candidate, size_t leg_id, double deltaR_Limit)
    {
        return FindMatchingTriggerObjects(descriptors, path_index, detail::GetTriggerObjectTypes(*candidate),
                                          candidate.GetMomentum(), leg_id, deltaR_Limit);
    }

    template<typename HiggsCandidate>
    void SetTriggerMatchBits(const analysis::TriggerDescriptors& descriptors,
                             const HiggsCandidate& candidate, double deltaR_Limit)
    {
        bool match_found = false;
        std::vector<analysis::TriggerDescriptors::PatternStruct> pattern_structs = descriptors.GetPatternStructs();
        
        for (size_t n = 0; n < pattern_structs.size(); ++n){
            analysis::TriggerDescriptors::PatternStruct pattern_struct = pattern_structs.at(n);
            
            //to be fixed
            std::map<size_t, TriggerObjectSet> matches;
            
            matches[first] = FindMatchingTriggerObjects(pattern_struct, n, candidate.GetFirstDaughter(),
                                                        deltaR_Limit);
            matches[second] = FindMatchingTriggerObjects(pattern_struct, n, candidate.GetSecondDaughter(),
                                                         deltaR_Limit);
            
            std::vector<const pat::TriggerObjectStandAlone*> comb_match;
            std::set_union(matches[1].begin(), matches[1].end(), matches[2].begin(), matches[2].end(),
                           std::back_inserter(comb_match));
            
            match_found = matches[1].size() >= 1 && matches[2].size() >= n_legs - 1 && comb_match.size() >= n_legs;
            
        }
        results.SetMatch(n, match_found); //n cos'è??

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
    analysis::TriggerDescriptors triggerDescriptors;
    std::vector<edm::ParameterSet> hltPaths;
    std::string pattern;
    std::vector<std::string> legs;
    std::map<CMSSW_Process, Handle<edm::TriggerResults>> triggerResultsMap;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;

    std::map<std::string,std::map<size_t,std::set<*pat::TriggerObjectStandAlone>>> path_legId_triggerObjPtr_map;
};

} // namespace analysis
