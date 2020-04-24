/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventCacheProvider.h"
#include "h-tautau/Analysis/include/TauUncertainties.h"
#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"


namespace analysis {
using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;
using LepCollection = std::vector<LepCandidate>;
using JetCandidate = Candidate<ntuple::TupleJet>;
using JetCollection = std::vector<JetCandidate>;
using FatJetCandidate = Candidate<ntuple::TupleFatJet>;
using FatJetCollection = std::vector<FatJetCandidate>;
using MET = MissingET<ntuple::TupleMet>;

class EventCandidate {
public:
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    EventCandidate(const ntuple::Event& _event, UncertaintySource _unc_source, UncertaintyScale _unc_scale);
    EventCandidate(const EventCandidate&) = delete;
    EventCandidate(EventCandidate&&) = delete;
    EventCandidate& operator=(const EventCandidate&) = delete;

    static void InitializeUncertainties(Period period, bool is_full, const std::string& working_path,
                                        TauIdDiscriminator tau_id_discriminator);

    static const jec::JECUncertaintiesWrapper& GetJecUncertainties();
    static const TauESUncertainties& GetTauESUncertainties();

    const LepCollection& GetLeptons() const;
    const JetCollection& GetJets() const;
    const FatJetCollection& GetFatJets() const;
    const MET& GetMET() const;
    bool IsSameAsCentral() const;
    const ntuple::Event& GetEvent() const;
    UncertaintySource GetUncSource() const;
    UncertaintyScale GetUncScale() const;
    UncertaintySource GetCacheUncSource() const;
    UncertaintyScale GetCacheUncScale() const;
    const EventIdentifier& GetEventId() const;
    Channel GetChannel() const;
    Period GetPeriod() const;

    const EventCacheProvider& GetCacheProvider();
    void SetCacheProvider(const std::shared_ptr<EventCacheProvider>& _cache_provider);
    void SetHHTagScores(size_t htt_index);

private:
    void CreateLeptons();
    void CreateJets();
    void CreateFatJets();

    Mutex mutex;
    const ntuple::Event* event;
    const UncertaintySource unc_source;
    const UncertaintyScale unc_scale;
    const EventIdentifier event_id;
    bool same_as_central;
    std::shared_ptr<EventCacheProvider> cache_provider;
    std::vector<ntuple::TupleLepton> tuple_leptons;
    std::vector<ntuple::TupleJet> tuple_jets;
    std::vector<ntuple::TupleFatJet> tuple_fatJets;
    ntuple::TupleMet tuple_met;
    std::vector<LepCandidate> lepton_candidates;
    std::vector<JetCandidate> jet_candidates;
    FatJetCollection fatJets;
    MET met;

    static const std::unique_ptr<boost::optional<jec::JECUncertaintiesWrapper>> jecUncertainties;
    static const std::unique_ptr<boost::optional<TauESUncertainties>> tauESUncertainties;
};

} // namespace analysis
