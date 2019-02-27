/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "AnalysisTools/Core/include/RootExt.h"

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"

#include "KinFitInterface.h"
#include "MT2.h"
#include "TriggerResults.h"

namespace analysis {

using ElectronCandidate = LeptonCandidate<ntuple::TupleElectron>;
using MuonCandidate = LeptonCandidate<ntuple::TupleMuon>;
using TauCandidate = LeptonCandidate<ntuple::TupleTau>;
using JetCandidate = Candidate<ntuple::TupleJet>;
using FatJetCandidate = Candidate<ntuple::TupleFatJet>;
using MET = MissingET<ntuple::TupleMet>;

namespace detail {
template<typename FirstLeg, typename SecondLeg>
constexpr Channel IdentifyChannel();

template<>
inline constexpr Channel IdentifyChannel<ElectronCandidate, TauCandidate>() { return Channel::ETau; }

template<>
inline constexpr Channel IdentifyChannel<MuonCandidate, TauCandidate>() { return Channel::MuTau; }

template<>
inline constexpr Channel IdentifyChannel<TauCandidate, TauCandidate>() { return Channel::TauTau; }

template<>
inline constexpr Channel IdentifyChannel<MuonCandidate, MuonCandidate>() { return Channel::MuMu; }
}

struct ChannelInfo {
    template<typename FirstLeg, typename SecondLeg>
    static constexpr Channel IdentifyChannel() { return detail::IdentifyChannel<FirstLeg, SecondLeg>(); }
};

template<int channel_id> struct ChannelLegInfo;
template<> struct ChannelLegInfo<static_cast<int>(Channel::ETau)> {
    using FirstLeg = ElectronCandidate; using SecondLeg = TauCandidate;
};
template<> struct ChannelLegInfo<static_cast<int>(Channel::MuTau)> {
    using FirstLeg = MuonCandidate; using SecondLeg = TauCandidate;
};
template<> struct ChannelLegInfo<static_cast<int>(Channel::TauTau)> {
    using FirstLeg = TauCandidate; using SecondLeg = TauCandidate;
};
template<> struct ChannelLegInfo<static_cast<int>(Channel::MuMu)> {
    using FirstLeg = MuonCandidate; using SecondLeg = MuonCandidate;
};

namespace jet_ordering {

    template<typename LorentzVector>
    struct JetInfo {
        LorentzVector p4;
        size_t index;
        double tag;

        JetInfo() : index(0), tag(0.0) { }

        JetInfo(const LorentzVector _p4, size_t _index, double _tag)
            : p4(_p4), index(_index), tag(_tag) { }

    };

    template<typename LorentzVector>
    bool CompareJets(const JetInfo<LorentzVector>& jet_1, const JetInfo<LorentzVector>& jet_2,
                     double pt_thr, double eta_thr)
    {
        const auto eta1 = std::abs(jet_1.p4.eta());
        const auto eta2 = std::abs(jet_2.p4.eta());
        if(eta1 < eta_thr && eta2 >= eta_thr) return true;
        if(eta1 >= eta_thr && eta2 < eta_thr) return false;
        const auto pt1 = jet_1.p4.pt();
        const auto pt2 = jet_2.p4.pt();
        if(pt1 > pt_thr && pt2 <= pt_thr) return true;
        if(pt1 <= pt_thr && pt2 > pt_thr) return false;

        if(jet_1.tag != jet_2.tag)
            return jet_1.tag > jet_2.tag;
        return pt1 > pt2;
    };

    template<typename LorentzVector>
    std::vector<JetInfo<LorentzVector>> OrderJets(
                                  const std::vector<JetInfo<LorentzVector>>& jet_info_vector,
                                  bool apply_hard_cut,
                                  double pt_cut = std::numeric_limits<double>::lowest(),
                                  double eta_cut = std::numeric_limits<double>::max())
    {
        const auto comparitor = [&](const JetInfo<LorentzVector>& jet_1,
                                    const JetInfo<LorentzVector>& jet_2) -> bool {
            return analysis::jet_ordering::CompareJets(jet_1, jet_2, pt_cut, eta_cut);
        };

        std::vector<JetInfo<LorentzVector>> jets_ordered;
        if(apply_hard_cut){
            for(size_t n = 0; n < jet_info_vector.size(); ++n) {
                if(jet_info_vector.at(n).p4.Pt() > pt_cut && std::abs(jet_info_vector.at(n).p4.eta()) < eta_cut)
                    jets_ordered.push_back(jet_info_vector.at(n));
            }
        }
        else
            jets_ordered = jet_info_vector;

        std::sort(jets_ordered.begin(), jets_ordered.end(), comparitor);
        return jets_ordered;
    }

}

class SummaryInfo {
public:
    using ProdSummary = ntuple::ProdSummary;

    explicit SummaryInfo(const ProdSummary& _summary, const std::string& _uncertainties_source = "");
    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors(Channel channel) const;
    const ProdSummary& operator*() const;
    const ProdSummary* operator->() const;
    const jec::JECUncertaintiesWrapper& GetJecUncertainties() const;

private:
    ProdSummary summary;
    std::map<Channel, std::shared_ptr<TriggerDescriptorCollection>> triggerDescriptors;
    std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;
};

class EventInfoBase {
public:
    using Event = ntuple::Event;
    using JetPair = ntuple::JetPair;
    using JetCollection = std::vector<JetCandidate>;
    using FatJetCollection = std::vector<FatJetCandidate>;
    using HiggsBBCandidate = CompositCandidate<JetCandidate, JetCandidate>;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    struct SelectedSignalJets{
        JetPair selectedBjetPair;
        JetPair selectedVBFjetPair;
        size_t n_bjets;

        SelectedSignalJets();
        bool HasBjetPair(size_t njets) const;
        bool HasVBFPair(size_t njets) const;
        bool isSelectedBjet(size_t n) const;
        bool isSelectedVBFjet(size_t n) const;
    };

    static SelectedSignalJets SelectSignalJets(const Event& event,
                                               const analysis::Period& period,
                                               JetOrdering jet_ordering);
    std::array<size_t,2> GetSelectedBjetIndices() const;
    std::set<size_t> GetSelectedBjetIndicesSet() const;

    EventInfoBase(const Event& _event, Period _period, JetOrdering _jet_ordering,
                  size_t _first_lepton_index, size_t _second_lepton_index,
                  const SummaryInfo* _summaryInfo = nullptr);

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

    virtual const AnalysisObject& GetLeg(size_t /*leg_id*/);
    virtual LorentzVector GetHiggsTTMomentum(bool /*useSVfit*/);

    size_t GetNJets() const;
    size_t GetNFatJets() const;
    const SelectedSignalJets& GetSelectedSignalJets() const;
    Period GetPeriod() const;
    JetOrdering GetJetOrdering() const;

    const JetCollection& GetJets();
    void SetJets(const JetCollection& new_jets);
    JetCollection SelectJets(double pt_cut = std::numeric_limits<double>::lowest(),
                             double eta_cut = std::numeric_limits<double>::max(),
                             JetOrdering jet_ordering = JetOrdering::DeepCSV,
                             const std::set<size_t>& jet_to_exclude_indexes = {});

    double GetHT(bool includeHbbJets, bool apply_pt_eta_cut);
    const FatJetCollection& GetFatJets();
    bool HasBjetPair() const;
    bool HasVBFjetPair() const;
    const JetCandidate& GetVBFJet(const size_t index);
    const JetCandidate& GetBJet(const size_t index);
    const HiggsBBCandidate& GetHiggsBB();
    const MET& GetMET();
    const size_t GetLegIndex(const size_t leg_id);

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

    LorentzVector GetResonanceMomentum(bool useSVfit, bool addMET);
    double GetMT2();
    const FatJetCandidate* SelectFatJet(double mass_cut, double deltaR_subjet_cut);
    void SetMvaScore(double _mva_score);
    double GetMvaScore() const;

    virtual std::shared_ptr<EventInfoBase> ApplyShiftBase(UncertaintySource uncertainty_source,
        UncertaintyScale scale) = 0;

protected:
    const Event* event;
    const SummaryInfo* summaryInfo;
    TriggerResults triggerResults;
    std::shared_ptr<Mutex> mutex;

private:
    template<typename LorentzVector>
    static bool PassEcalNoiceVetoJets(const LorentzVector& jet_p4, Period period)
    {
        if(period !=  analysis::Period::Run2017)
            return true;

        const double abs_eta = std::abs(jet_p4.eta());
        return !(jet_p4.pt() < cuts::hh_bbtautau_2017::jetID::max_pt_veto &&
                    abs_eta > cuts::hh_bbtautau_2017::jetID::eta_low_veto &&
                    abs_eta < cuts::hh_bbtautau_2017::jetID::eta_high_veto);
    }

private:
    EventIdentifier eventIdentifier;
    SelectedSignalJets selected_signal_jets;
    Period period;
    JetOrdering jet_ordering;
    size_t first_lepton_index;
    size_t second_lepton_index;

    std::shared_ptr<std::list<ntuple::TupleJet>> tuple_jets;
    std::shared_ptr<JetCollection> jets;
    std::shared_ptr<std::list<ntuple::TupleFatJet>> tuple_fatJets;
    std::shared_ptr<FatJetCollection> fatJets;
    std::shared_ptr<HiggsBBCandidate> higgs_bb;
    std::shared_ptr<ntuple::TupleMet> tuple_met;
    std::shared_ptr<MET> met;
    std::shared_ptr<kin_fit::FitResults> kinfit_results;
    boost::optional<double> mt2;
    double mva_score;

};

template<typename _FirstLeg, typename _SecondLeg>
class EventInfo : public EventInfoBase {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using FirstTupleLeg = typename FirstLeg::PATObject;
    using SecondTupleLeg = typename SecondLeg::PATObject;
    using HiggsTTCandidate = CompositCandidate<FirstLeg, SecondLeg>;

    static constexpr Channel channel = ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>();

    EventInfo(const Event& _event, Period _period, JetOrdering _jet_ordering,
              size_t _first_lepton_index, size_t _second_lepton_index,
              const SummaryInfo* _summaryInfo = nullptr) :
        EventInfoBase(_event, _period, _jet_ordering, _first_lepton_index,
            _second_lepton_index, _summaryInfo)
    {
        if(summaryInfo)
            triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors(channel));
    }

    using EventInfoBase::EventInfoBase;

    const FirstLeg& GetFirstLeg()
    {
        Lock lock(*mutex);
        if(!leg1) {
            tuple_leg1 = std::shared_ptr<FirstTupleLeg>(new FirstTupleLeg(*event, GetLegIndex(1)));
            leg1 = std::shared_ptr<FirstLeg>(new FirstLeg(*tuple_leg1, tuple_leg1->iso()));
        }
        return *leg1;
    }

    const SecondLeg& GetSecondLeg()
    {
        Lock lock(*mutex);
        if(!leg2) {
            tuple_leg2 = std::shared_ptr<SecondTupleLeg>(new SecondTupleLeg(*event, GetLegIndex(2)));
            leg2 = std::shared_ptr<SecondLeg>(new SecondLeg(*tuple_leg2, tuple_leg2->iso()));
        }
        return *leg2;
    }

    virtual const AnalysisObject& GetLeg(size_t leg_id) override
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
                higgs_tt_sv = std::shared_ptr<HiggsTTCandidate>(
                              new HiggsTTCandidate(GetFirstLeg(), GetSecondLeg(), event->SVfit_p4));
            }
            return *higgs_tt_sv;
        }
        if(!higgs_tt)
            higgs_tt = std::shared_ptr<HiggsTTCandidate>(new HiggsTTCandidate(GetFirstLeg(), GetSecondLeg()));
        return *higgs_tt;
    }

    virtual LorentzVector GetHiggsTTMomentum(bool useSVfit) override
    {
        return GetHiggsTT(useSVfit).GetMomentum();
    }

    virtual std::shared_ptr<EventInfoBase> ApplyShiftBase(UncertaintySource uncertainty_source,
        UncertaintyScale scale) override
    {
        EventInfo<FirstLeg, SecondLeg> event_info = ApplyShift(uncertainty_source, scale);
        return std::make_shared<EventInfo<FirstLeg, SecondLeg>>(std::move(event_info));
    }

    EventInfo<FirstLeg, SecondLeg> ApplyShift(UncertaintySource uncertainty_source,
        UncertaintyScale scale)
    {
        EventInfo<FirstLeg, SecondLeg> shifted_event_info(*this);
        const SummaryInfo& summaryInfo = shifted_event_info.GetSummaryInfo();
        const jec::JECUncertaintiesWrapper& jecUncertainties = summaryInfo.GetJecUncertainties();
        const JetCollection& jets = shifted_event_info.GetJets();
        const auto& other_jets_p4 = event->other_jets_p4;
        auto shifted_met_p4(shifted_event_info.GetMET().GetMomentum());
        const JetCollection& corrected_jets = jecUncertainties.ApplyShift(jets,uncertainty_source,scale,&other_jets_p4,&shifted_met_p4);
        shifted_event_info.SetJets(corrected_jets);
        shifted_event_info.SetMetMomentum(shifted_met_p4);
        return shifted_event_info;
    }

private:
    std::shared_ptr<FirstTupleLeg> tuple_leg1;
    std::shared_ptr<FirstLeg> leg1;
    std::shared_ptr<SecondTupleLeg> tuple_leg2;
    std::shared_ptr<SecondLeg> leg2;
    std::shared_ptr<HiggsTTCandidate> higgs_tt, higgs_tt_sv;
};

std::shared_ptr<EventInfoBase> MakeEventInfo(Channel channel, const EventInfoBase::Event& event,
                                             Period period, JetOrdering jet_ordering,
                                             const SummaryInfo* summaryInfo = nullptr);
} // namespace analysis
