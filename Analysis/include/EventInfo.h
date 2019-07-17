/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "AnalysisTools/Core/include/RootExt.h"

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TriggerResults.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "h-tautau/Cuts/include/H_tautau_2017_baseline.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"

#include "SVfitAnaInterface.h"
#include "KinFitInterface.h"
#include "MT2.h"

#include "SignalObjectSelector.h"

#include <numeric>


namespace analysis {

using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;
using JetCandidate = Candidate<ntuple::TupleJet>;
using FatJetCandidate = Candidate<ntuple::TupleFatJet>;
using MET = MissingET<ntuple::TupleMet>;



class SummaryInfo {
public:
    using ProdSummary = ntuple::ProdSummary;

    explicit SummaryInfo(const ProdSummary& _summary, const Channel& _channel, const std::string& _uncertainties_source = "",
                         const std::string& _trigger_cfg = "");
    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors() const;
    const ProdSummary& operator*() const;
    const ProdSummary* operator->() const;
    const jec::JECUncertaintiesWrapper& GetJecUncertainties() const;

private:
    ProdSummary summary;
    std::shared_ptr<const TriggerDescriptorCollection> triggerDescriptors;
    std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;

};

class EventInfoBase {
public:
    using Event = ntuple::Event;
    using JetPair = ntuple::JetPair;
    using JetCollection = std::vector<JetCandidate>;
    using FatJetCollection = std::vector<FatJetCandidate>;
    using HiggsBBCandidate = CompositeCandidate<JetCandidate, JetCandidate>;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;
    using HiggsTTCandidate = CompositeCandidate<LepCandidate, LepCandidate>;


    std::array<size_t,2> GetSelectedBjetIndices() const;
    std::set<size_t> GetSelectedBjetIndicesSet() const;

    Channel GetChannel() const { return static_cast<Channel>(event->channelId); }

    EventInfoBase(const Event& _event, const SummaryInfo* _summaryInfo,
                  size_t _selected_htt_index, const SignalObjectSelector::SelectedSignalJets& _selected_signal_jets,
                  Period _period, JetOrdering _jet_ordering);


    EventInfoBase(const EventInfoBase& ) = default; //copy constructor
    virtual ~EventInfoBase(){} //destructor

    EventInfoBase& operator= ( const EventInfoBase& ) = default; //assignment


    const Event& operator*() const;
    const Event* operator->() const;

    const EventIdentifier& GetEventId() const;
    EventEnergyScale GetEnergyScale() const;
    const TriggerResults& GetTriggerResults() const;
    const SummaryInfo& GetSummaryInfo() const;
    static const kin_fit::FitProducer& GetKinFitProducer();
    static const sv_fit_ana::FitProducer& GetSVFitProducer();

    // virtual const Candidate& GetLeg(size_t /*leg_id*/);
    // virtual LorentzVector GetHiggsTTMomentum(bool /*useSVfit*/);

    size_t GetNJets() const;
    size_t GetNFatJets() const;
    const SignalObjectSelector::SelectedSignalJets& GetSelectedSignalJets() const;
    Period GetPeriod() const;
    JetOrdering GetJetOrdering() const;

    const JetCollection& GetJets();
    void SetJets(const JetCollection& new_jets);
    JetCollection SelectJets(double pt_cut = std::numeric_limits<double>::lowest(),
                             double eta_cut = std::numeric_limits<double>::max(),
                             bool applyPu = false, bool passBtag = false,
                             JetOrdering jet_ordering = JetOrdering::DeepCSV,
                             const std::set<size_t>& jet_to_exclude_indexes = {},
                             double low_eta_cut = 0);

    double GetHT(bool includeHbbJets, bool apply_pt_eta_cut);
    const FatJetCollection& GetFatJets();
    bool HasBjetPair() const;
    bool HasVBFjetPair() const;
    const JetCandidate& GetVBFJet(const size_t index);
    const JetCandidate& GetBJet(const size_t index);
    const HiggsBBCandidate& GetHiggsBB();
    const MET& GetMET();
    size_t GetLegIndex(const size_t leg_id);
    static bool PassDefaultLegSelection(const ntuple::TupleLepton& lepton, Channel channel);

    template<typename LorentzVector>
    void SetMetMomentum(const LorentzVector& new_met_p4)
    {
        Lock lock(*mutex);
        if(!tuple_met)
            tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
        met = std::make_shared<MET>(*tuple_met, tuple_met->cov());
        met->SetMomentum(new_met_p4);
    }

    const kin_fit::FitResults& GetKinFitResults();
    const sv_fit_ana::FitResults& GetSVFitResults();

    LorentzVector GetResonanceMomentum(bool useSVfit, bool addMET);
    double GetMT2();
    const FatJetCandidate* SelectFatJet(double mass_cut, double deltaR_subjet_cut);
    void SetMvaScore(double _mva_score);
    double GetMvaScore() const;

    const LepCandidate& GetFirstLeg();
    const LepCandidate& GetSecondLeg();
    std::shared_ptr<EventInfoBase> ApplyShift(UncertaintySource uncertainty_source, UncertaintyScale scale);

    const LepCandidate& GetLeg(size_t leg_id)
    {
        if(leg_id == 1) return GetFirstLeg();
        if(leg_id == 2) return GetSecondLeg();
        throw exception("Invalid leg id = %1%.") % leg_id;
    }

    const HiggsTTCandidate& GetHiggsTT(bool useSVfit)
    {
        Lock lock(*mutex);
        if(useSVfit) {
            if(!higgs_tt_sv) {
                if(!GetSVFitResults().has_valid_momentum) throw exception("SVFit not converged");
                higgs_tt_sv = std::make_shared<HiggsTTCandidate>(GetFirstLeg(), GetSecondLeg(), GetSVFitResults().momentum);
            }
            return *higgs_tt_sv;
        }
        if(!higgs_tt)
            higgs_tt = std::make_shared<HiggsTTCandidate>(GetFirstLeg(), GetSecondLeg());
        return *higgs_tt;
    }

    LorentzVector GetHiggsTTMomentum(bool useSVfit)
    {
        return GetHiggsTT(useSVfit).GetMomentum();
    }






protected:
    const Event* event;
    const SummaryInfo* summaryInfo;
    TriggerResults triggerResults;
    std::shared_ptr<Mutex> mutex;
    size_t selected_htt_index;

private:
    EventIdentifier eventIdentifier;
    SignalObjectSelector::SelectedSignalJets selected_signal_jets;
    Period period;
    JetOrdering jet_ordering;


    std::shared_ptr<std::list<ntuple::TupleJet>> tuple_jets;
    std::shared_ptr<JetCollection> jets;
    std::shared_ptr<std::list<ntuple::TupleFatJet>> tuple_fatJets;
    std::shared_ptr<FatJetCollection> fatJets;
    std::shared_ptr<HiggsBBCandidate> higgs_bb;
    std::shared_ptr<ntuple::TupleMet> tuple_met;
    std::shared_ptr<MET> met;
    std::shared_ptr<kin_fit::FitResults> kinfit_results;
    std::shared_ptr<sv_fit_ana::FitResults> svfit_results;
    boost::optional<double> mt2;
    double mva_score;
    std::shared_ptr<ntuple::TupleLepton> tuple_leg1;
    std::shared_ptr<LepCandidate> leg1;
    std::shared_ptr<ntuple::TupleLepton> tuple_leg2;
    std::shared_ptr<LepCandidate> leg2;
    std::shared_ptr<HiggsTTCandidate> higgs_tt, higgs_tt_sv;

};

boost::optional<EventInfoBase> CreateEventInfo(const ntuple::Event& event,
                                               const SignalObjectSelector& signalObjectSelector,
                                               const SummaryInfo* summaryInfo = nullptr,
                                               Period period = analysis::Period::Run2017,
                                               JetOrdering jet_ordering = JetOrdering::DeepCSV);

} // namespace analysis
