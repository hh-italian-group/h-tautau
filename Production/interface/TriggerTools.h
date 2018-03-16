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
#include "h-tautau/Production/interface/TriggerFileDescriptor.h"
#include "h-tautau/Production/interface/TriggerFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"

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
LegType GetTriggerObjectTypes<pat::Electron>(const pat::Electron&) { return LegType::e; }

template<>
LegType GetTriggerObjectTypes<pat::Muon>(const pat::Muon&) { return LegType::mu; }

template<>
LegType GetTriggerObjectTypes<pat::Tau>(const pat::Tau&) { return LegType::tau; }

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
                 std::string triggerCfg, Channel channel);

    static trigger_tools::TriggerFileDescriptorCollection ReadConfig(const std::string& cfg_path,
                                                                    trigger_tools::SetupDescriptor& setup);

    static TriggerDescriptorCollection CreateTriggerDescriptors(const trigger_tools::TriggerFileDescriptorCollection& trigger_file_descriptors, Channel channel);

    void Initialize(const edm::Event& iEvent);

    void SetTriggerAcceptBits(TriggerResults& results);

    std::vector<TriggerObjectSet> FindMatchingTriggerObjects(size_t desc_index,
            const LorentzVector& candidateMomentum, LegType candidate_type, double deltaR_Limit);


    template<typename HiggsCandidate>
    void SetTriggerMatchBits(TriggerResults& results, const HiggsCandidate& candidate, double deltaR_Limit)
    { 
        std::array<std::vector<TriggerTools::TriggerObjectSet>, 2> array_matched_legId_triggerObjectSet;
        for (size_t n = 0; n < triggerDescriptors.size(); ++n){
            const auto& descriptor = triggerDescriptors.at(n);
            std::vector<TriggerObjectSet> matches_first, matches_second;
            matches_first = FindMatchingTriggerObjects(n,candidate.GetFirstDaughter().GetMomentum(),
                                                        detail::GetTriggerObjectTypes(*candidate.GetFirstDaughter()),
                                                        deltaR_Limit);
            matches_second = FindMatchingTriggerObjects(n,candidate.GetSecondDaughter().GetMomentum(),
                                                         detail::GetTriggerObjectTypes(*candidate.GetSecondDaughter()),
                                                         deltaR_Limit);
            array_matched_legId_triggerObjectSet.at(0) = matches_first;
            array_matched_legId_triggerObjectSet.at(1) = matches_second;
            results.SetMatch(n, TriggerMatchFound(array_matched_legId_triggerObjectSet, descriptor.legs_info.size()));
        }
    }

    bool TriggerMatchFound(const std::array<std::vector<TriggerTools::TriggerObjectSet>, 2>& array_matched_legId_triggerObjectSet,
                           const size_t n_legs_total)
    {
        if(n_legs_total == 0) return true;

        bool match_found = false;
        for(size_t flip = 0; !match_found && flip < array_matched_legId_triggerObjectSet.size(); ++flip) {
            const size_t first = (flip % 2) + 1, second = ((flip + 1) % 2) + 1;
            if(n_legs_total == 1)
                return array_matched_legId_triggerObjectSet.at(flip).at(first).size() >= n_legs_total ||
                        array_matched_legId_triggerObjectSet.at(flip).at(second).size() >= n_legs_total;
            if(n_legs_total == 2){
                std::vector<const pat::TriggerObjectStandAlone*> comb_match;
                std::set_union(array_matched_legId_triggerObjectSet.at(flip).at(first).begin(),
                               array_matched_legId_triggerObjectSet.at(flip).at(first).end(),
                               array_matched_legId_triggerObjectSet.at(flip).at(second).begin(),
                               array_matched_legId_triggerObjectSet.at(flip).at(second).end(),
                                std::back_inserter(comb_match));

                match_found = array_matched_legId_triggerObjectSet.at(flip).at(0).size() >= 1 &&
                        array_matched_legId_triggerObjectSet.at(flip).at(1).size() >= n_legs_total - 1 &&
                        comb_match.size() >= n_legs_total;
                if(match_found) return true;
            }
       }
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
    std::map<LegType, double> deltaPt_map;
    std::map<CMSSW_Process, Handle<edm::TriggerResults>> triggerResultsMap;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;

    std::vector<std::vector<TriggerObjectSet>> path_legId_triggerObjPtr_vector;
};

} // namespace analysis
