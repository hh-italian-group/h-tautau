/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TriggerResults.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Analysis/include/EventCandidate.h"
#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "h-tautau/Analysis/include/MT2.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "h-tautau/Analysis/include/SVfitAnaInterface.h"

namespace analysis {

class SummaryInfo {
public:
    using ProdSummary = ntuple::ProdSummary;

    explicit SummaryInfo(const ProdSummary& _summary, const Channel& _channel,
                         const std::string& _trigger_cfg = "", const std::vector<std::string>& _triggers = {},
                         const std::vector<std::string>& _vbf_triggers = {});
    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors() const;
    const ProdSummary& operator*() const;
    const ProdSummary* operator->() const;
    const jec::JECUncertaintiesWrapper& GetJecUncertainties() const;
    const std::vector<std::string>& GetTriggers() const;
    const std::vector<std::string>& GetVbfTriggers() const;

private:
    ProdSummary summary;
    std::shared_ptr<const TriggerDescriptorCollection> triggerDescriptors;
    std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;
    std::vector<std::string> triggers, vbf_triggers;
};

class EventInfo {
public:
    using Event = ntuple::Event;
    using LegPair = ntuple::LegPair;
    using SummaryInfoPtr = std::shared_ptr<const SummaryInfo>;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;
    using HiggsTTCandidate = CompositeCandidate<LepCandidate, LepCandidate>;
    using HiggsBBCandidate = CompositeCandidate<JetCandidate, JetCandidate>;
    using SelectedSignalJets = SignalObjectSelector::SelectedSignalJets;

    EventInfo(const std::shared_ptr<EventCandidate>& _event_candidate,
              const std::shared_ptr<const SummaryInfo>& _summaryInfo, size_t _selected_htt_index,
              const SelectedSignalJets& _selected_signal_jets, const BTagger& bTagger);

    EventInfo(const EventInfo&) = delete;
    EventInfo& operator=(const EventInfo&) = delete;

    const Event& operator*() const;
    const Event* operator->() const;

    const EventCandidate& GetEventCandidate() const;
    const SummaryInfo& GetSummaryInfo() const;
    size_t GetHttIndex() const;
    const ntuple::LegPair& GetSelectedTauIndices() const;
    const SelectedSignalJets& GetSelectedSignalJets() const;
    const BTagger& GetBTagger() const;
    const TriggerResults& GetTriggerResults() const;

    // EventCandidate proxies
    Channel GetChannel() const;
    Period GetPeriod() const;
    const EventIdentifier& GetEventId() const;
    const MET& GetMET() const;

    const LepCandidate& GetLeg(size_t leg_id) const;
    const LepCandidate& GetFirstLeg() const;
    const LepCandidate& GetSecondLeg() const;
    const HiggsTTCandidate& GetHiggsTT(bool useSVfit, bool allow_calc = false);
    boost::optional<LorentzVector> GetHiggsTTMomentum(bool useSVfit, bool allow_calc = false);

    bool HasBjetPair() const;
    const JetCandidate& GetBJet(size_t index) const;
    const HiggsBBCandidate& GetHiggsBB();

    bool HasVBFjetPair() const;
    const JetCandidate& GetVBFJet(size_t index) const;

    boost::optional<LorentzVector> GetResonanceMomentum(bool useSVfit, bool addMET, bool allow_calc = false);
    double GetHT(bool includeHbbJets, bool includeVBFJets) const;
    const sv_fit_ana::FitResults& GetSVFitResults(bool allow_calc = false, int verbosity = 0);
    const kin_fit::FitResults& GetKinFitResults(bool allow_calc = false, int verbosity = 0);
    double GetMT2();

    const std::vector<const JetCandidate*>& GetCentralJets();
    const std::vector<const JetCandidate*>& GetForwardJets();
    const std::vector<const JetCandidate*>& GetAllJets();

    void SetMvaScore(double _mva_score);
    double GetMvaScore() const;

    bool PassNormalTriggers();
    bool PassVbfTriggers();

    [[ noreturn ]] void ThrowException(const std::string& message) const;

    static std::unique_ptr<EventInfo> Create(const ntuple::Event& event,
                                             const SignalObjectSelector& signalObjectSelector,
                                             const BTagger& bTagger, DiscriminatorWP btag_wp,
                                             SummaryInfoPtr summaryInfo = SummaryInfoPtr(),
                                             UncertaintySource unc_source = UncertaintySource::None,
                                             UncertaintyScale unc_scale = UncertaintyScale::Central,
                                             bool is_sync = false, bool debug = false);

    static std::unique_ptr<EventInfo> Create(const std::shared_ptr<EventCandidate>& event_candidate,
                                             const SignalObjectSelector& signalObjectSelector,
                                             const BTagger& bTagger, DiscriminatorWP btag_wp,
                                             SummaryInfoPtr summaryInfo = SummaryInfoPtr(),
                                             bool is_sync = false, bool debug = false);

private:
    std::vector<std::string> FilterTriggers(const std::vector<std::string>& trigger_names) const;

private:
    Mutex mutex;
    const std::shared_ptr<EventCandidate> event_candidate;
    const SummaryInfoPtr summaryInfo;
    const TriggerResults triggerResults;
    const BTagger bTagger;
    const size_t selected_htt_index;
    const ntuple::LegPair selected_tau_indices;
    const SelectedSignalJets selected_signal_jets;

    boost::optional<HiggsTTCandidate> higgs_tt, higgs_tt_sv;
    boost::optional<HiggsBBCandidate> higgs_bb;
    boost::optional<kin_fit::FitResults> kinfit_results;
    boost::optional<sv_fit_ana::FitResults> svfit_results;
    boost::optional<double> mt2, mva_score;
    boost::optional<std::vector<const JetCandidate*>> central_jets, forward_jets, all_jets;
    boost::optional<bool> pass_triggers, pass_vbf_triggers;
};

} // namespace analysis
