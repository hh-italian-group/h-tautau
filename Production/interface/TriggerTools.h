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
inline LegType GetTriggerObjectTypes(const PatObject&);

template<>
inline LegType GetTriggerObjectTypes<pat::Electron>(const pat::Electron&) { return LegType::e; }

template<>
inline LegType GetTriggerObjectTypes<pat::Muon>(const pat::Muon&) { return LegType::mu; }

template<>
inline LegType GetTriggerObjectTypes<pat::Tau>(const pat::Tau&) { return LegType::tau; }

} // namespace detail

class TriggerTools {
public:
    using L1ParticlePtrSet = std::set<const l1extra::L1JetParticle*>;
    template<typename T> using EDGetTokenT = edm::EDGetTokenT<T>;
    template<typename T> using Handle = edm::Handle<T>;
    using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;
    using VectorTriggerObjectSet = std::vector<TriggerTools::TriggerObjectSet>;

    TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                 EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                 EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                 EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token,
                 const std::string& triggerCfg, Channel channel);

    static trigger_tools::TriggerFileDescriptorCollection ReadConfig(const std::string& cfg_path,
                                                                    trigger_tools::SetupDescriptor& setup);

    static TriggerDescriptorCollection CreateTriggerDescriptors(const trigger_tools::TriggerFileDescriptorCollection& trigger_file_descriptors, Channel channel);

    void Initialize(const edm::Event& iEvent);

    void SetTriggerAcceptBits(TriggerResults& results);

    VectorTriggerObjectSet FindMatchingTriggerObjects(size_t desc_index,
            const LorentzVector& candidateMomentum, LegType candidate_type, double deltaR_Limit) const;


    template<typename HiggsCandidate>
    void SetTriggerMatchBits(TriggerResults& results, const HiggsCandidate& candidate, double deltaR_Limit)
    { 
        std::array<VectorTriggerObjectSet, 2> matched_legIds;
        for (size_t n = 0; n < triggerDescriptors.size(); ++n){
            const auto& descriptor = triggerDescriptors.at(n);
            matched_legIds.at(0) = FindMatchingTriggerObjects(n,candidate.GetFirstDaughter().GetMomentum(),
                                                              detail::GetTriggerObjectTypes(*candidate.GetFirstDaughter()),
                                                              deltaR_Limit);
            matched_legIds.at(1) = FindMatchingTriggerObjects(n,candidate.GetSecondDaughter().GetMomentum(),
                                                              detail::GetTriggerObjectTypes(*candidate.GetSecondDaughter()),
                                                              deltaR_Limit);
            results.SetMatch(n, TriggerMatchFound(matched_legIds, descriptor.legs_info.size()));
        }
    }

    bool TriggerMatchFound(const std::array<VectorTriggerObjectSet, 2>& matched_legIds,
                           const size_t n_legs_total);


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

    std::vector<VectorTriggerObjectSet> pathTriggerObjects;
};

} // namespace analysis
