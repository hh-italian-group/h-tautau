/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <list>
#include <algorithm>
#include "EventTuple.h"
#include "AnalysisTypes.h"
#include "RootExt.h"
#include "KinFitInterface.h"
#include "Candidate.h"
#include "TupleObjects.h"
#include "TriggerResults.h"
#include "SummaryTuple.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "hh-bbtautau/Analysis/include/MT2.h"

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


enum class JetOrdering { NoOrdering, Pt, CSV, DeepCSV };

class SummaryInfo {
public:
    using ProdSummary = ntuple::ProdSummary;

    explicit SummaryInfo(const ProdSummary& _summary) : summary(_summary)
    {
        for(size_t n = 0; n < summary.triggers_channel.size(); ++n) {
            const int channel_id = summary.triggers_channel.at(n);
            const Channel channel = static_cast<Channel>(channel_id);
            if(!triggerDescriptors.count(channel))
                triggerDescriptors[channel] = std::make_shared<TriggerDescriptorCollection>();
            triggerDescriptors[channel]->Add(summary.triggers_pattern.at(n), {});
        }
    }

    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors(Channel channel) const
    {
        if(!triggerDescriptors.count(channel))
            throw exception("Information for channel %1% not found.") % channel;
        return triggerDescriptors.at(channel);
    }

    const ProdSummary& operator*() const { return summary; }
    const ProdSummary* operator->() const { return &summary; }

private:
    ProdSummary summary;
    std::map<Channel, std::shared_ptr<TriggerDescriptorCollection>> triggerDescriptors;
};

class EventInfoBase {
public:
    using Event = ntuple::Event;
    using JetPair = ntuple::JetPair;
    using JetCollection = std::vector<JetCandidate>;
    using FatJetCollection = std::vector<FatJetCandidate>;
    using HiggsBBCandidate = CompositCandidate<JetCandidate, JetCandidate>;

    static JetPair SelectBjetPair(const Event& event, double pt_cut = std::numeric_limits<double>::lowest(),
                                   double eta_cut = std::numeric_limits<double>::max(),
                                   JetOrdering jet_ordering = JetOrdering::CSV)
    {
        const auto orderer = [&](size_t j1, size_t j2) -> bool {
            if(jet_ordering == JetOrdering::Pt)
                return event.jets_p4.at(j1).Pt() > event.jets_p4.at(j2).Pt();
            if(jet_ordering == JetOrdering::CSV){
                if(event.jets_csv.at(j1) != event.jets_csv.at(j2))
                    return event.jets_csv.at(j1) > event.jets_csv.at(j2);
                return event.jets_p4.at(j1).Pt() > event.jets_p4.at(j2).Pt();
            }
            if(jet_ordering == JetOrdering::DeepCSV){
                if(event.jets_deepCsv_b.at(j1)+event.jets_deepCsv_bb.at(j1) != event.jets_deepCsv_b.at(j2)+event.jets_deepCsv_bb.at(j2))
                    return event.jets_deepCsv_b.at(j1)+event.jets_deepCsv_bb.at(j1) > event.jets_deepCsv_b.at(j2)+event.jets_deepCsv_bb.at(j2);
                return event.jets_p4.at(j1).Pt() > event.jets_p4.at(j2).Pt();
            }

            throw exception("Unsupported jet ordering for b-jet pair selection.");
        };

        std::vector<size_t> indexes;
        for(size_t n = 0; n < event.jets_p4.size(); ++n) {
            if(event.jets_p4.at(n).Pt() > pt_cut && std::abs(event.jets_p4.at(n).eta()) < eta_cut)
                indexes.push_back(n);
        }
        std::sort(indexes.begin(), indexes.end(), orderer);
        JetPair selected_pair = ntuple::UndefinedJetPair();
        if(indexes.size() >= 1)
            selected_pair.first = indexes.at(0);
        if(indexes.size() >= 2)
            selected_pair.second = indexes.at(1);
        return selected_pair;
    }

    std::array<size_t,2> GetSelectedBjetIndices() const
    {
        std::array<size_t,2> bjet_indexes;
        bjet_indexes[0] = selected_bjet_pair.first;
        bjet_indexes[1] = selected_bjet_pair.second;
        return bjet_indexes;
    }

    std::set<size_t> GetSelectedBjetIndicesSet() const
    {
        std::set<size_t> bjet_indexes;
        bjet_indexes.insert(selected_bjet_pair.first);
        bjet_indexes.insert(selected_bjet_pair.second);
        return bjet_indexes;
    }


    static constexpr int verbosity = 0;

    EventInfoBase(const Event& _event, const JetPair& _selected_bjet_pair = JetPair(0, 1),
                  const SummaryInfo* _summaryInfo = nullptr) :
        event(&_event), summaryInfo(_summaryInfo), eventIdentifier(_event.run, _event.lumi, _event.evt),
        selected_bjet_pair(_selected_bjet_pair),
        has_bjet_pair(selected_bjet_pair.first < GetNJets() && selected_bjet_pair.second < GetNJets())
    {
        triggerResults.SetAcceptBits(event->trigger_accepts);
        triggerResults.SetMatchBits(event->trigger_matches);
    }

    virtual ~EventInfoBase() {}

    const Event& operator*() const { return *event; }
    const Event* operator->() const { return event; }

    const EventIdentifier& GetEventId() const { return eventIdentifier; }
    EventEnergyScale GetEnergyScale() const { return static_cast<EventEnergyScale>(event->eventEnergyScale); }
    const TriggerResults& GetTriggerResults() const { return triggerResults; }
    const SummaryInfo& GetSummaryInfo() const
    {
        if(!summaryInfo)
            throw exception("SummaryInfo was not provided for this event.");
        return *summaryInfo;
    }

    virtual const AnalysisObject& GetLeg(size_t /*leg_id*/) { throw exception("Method not supported."); }
    virtual LorentzVector GetHiggsTTMomentum(bool /*useSVfit*/) { throw exception("Method not supported."); }

    size_t GetNJets() const { return event->jets_p4.size(); }
    size_t GetNFatJets() const { return event->fatJets_p4.size(); }

    const JetCollection& GetJets()
    {
        if(!jets) {
            jets = std::shared_ptr<JetCollection>(new JetCollection());
            for(size_t n = 0; n < GetNJets(); ++n) {
                tuple_jets.push_back(ntuple::TupleJet(*event, n));
                jets->push_back(JetCandidate(tuple_jets.back()));
            }
        }
        return *jets;
    }

    JetCollection SelectJets(double pt_cut = std::numeric_limits<double>::lowest(),
                             double eta_cut = std::numeric_limits<double>::max(),
                             double csv_cut = std::numeric_limits<double>::lowest(),
                             JetOrdering jet_ordering = JetOrdering::CSV,
                             const std::set<size_t>& bjet_indexes = {})
    {
        const auto orderer = [&](const JetCandidate& j1, const JetCandidate& j2) -> bool {
            if(jet_ordering == JetOrdering::Pt)
                return j1.GetMomentum().Pt() > j2.GetMomentum().Pt();
            if(jet_ordering == JetOrdering::DeepCSV){
                if(j1->deepcsv() != j2->deepcsv())
                    return j1->deepcsv() > j2->deepcsv();
                return j1.GetMomentum().Pt() > j2.GetMomentum().Pt();
            }
            if(jet_ordering == JetOrdering::CSV){
                if(j1->csv() != j2->csv())
                    return j1->csv() > j2->csv();
                return j1.GetMomentum().Pt() > j2.GetMomentum().Pt();
            }
            throw exception("Unsupported jet ordering for jet selection.");
        };

        const JetCollection& all_jets = GetJets();
        JetCollection selected_jets;
        for(size_t n = 0; n < all_jets.size(); ++n) {
            if(bjet_indexes.count(n)) continue;
            const JetCandidate& jet = all_jets.at(n);
            if(jet.GetMomentum().Pt() > pt_cut && std::abs(jet.GetMomentum().eta()) < eta_cut
                    && jet->csv() > csv_cut)
                selected_jets.push_back(jet);
        }

        if(jet_ordering != JetOrdering::NoOrdering)
            std::sort(selected_jets.begin(), selected_jets.end(), orderer);
        return selected_jets;
    }

    static JetPair SelectVBFJetPair(const JetCollection& jets_vbf)
    {
        double max_mjj = -std::numeric_limits<double>::infinity();
        JetPair selected_pair = ntuple::UndefinedJetPair();
        for(size_t n = 0; n < jets_vbf.size(); ++n) {
            for(size_t h = n+1; h < jets_vbf.size(); ++h) {
                const JetCandidate& jet_1 = jets_vbf.at(n);
                const JetCandidate& jet_2 = jets_vbf.at(h);
                const LorentzVector jet_12 = jet_1.GetMomentum() + jet_2.GetMomentum();
                if(jet_12.M() > max_mjj){
                    max_mjj = jet_12.M();
                    if(jet_1.GetMomentum().Pt() > jet_2.GetMomentum().Pt())
                        selected_pair = std::make_pair(n,h);
                    else
                        selected_pair = std::make_pair(h,n);
                }
            }
        }

        return selected_pair;
    }

    double GetHT(bool includeHbbJets = true)
    {
        const JetCollection& all_jets = GetJets();
        const std::set<size_t> bjet_indexes = GetSelectedBjetIndicesSet();
        double sum = 0;
        for(unsigned n = 0; n < all_jets.size(); ++n) {
            if(!includeHbbJets && bjet_indexes.count(n)) continue;
            const JetCandidate& jet = all_jets.at(n);
            sum += jet.GetMomentum().Pt();
        }
        return sum;
    }


    const FatJetCollection& GetFatJets()
    {
        if(!fatJets) {
            fatJets = std::shared_ptr<FatJetCollection>(new FatJetCollection());
            for(size_t n = 0; n < GetNFatJets(); ++n) {
                tuple_fatJets.push_back(ntuple::TupleFatJet(*event, n));
                fatJets->push_back(FatJetCandidate(tuple_fatJets.back()));
            }
        }
        return *fatJets;
    }

    bool HasBjetPair() const { return has_bjet_pair; }

    const HiggsBBCandidate& GetHiggsBB()
    {
        if(!HasBjetPair())
            throw exception("Can't create H->bb candidate.");
        if(!higgs_bb) {
            const auto& jets = GetJets();
            higgs_bb = std::shared_ptr<HiggsBBCandidate>(
                       new HiggsBBCandidate(jets.at(selected_bjet_pair.first), jets.at(selected_bjet_pair.second)));
        }
        return *higgs_bb;
    }

    const MET& GetMET()
    {
        if(!met) {
            tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
            met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
        }
        return *met;
    }

    const kin_fit::FitResults& GetKinFitResults()
    {
        if(!HasBjetPair())
            throw exception("Can't retrieve KinFit results.");
        if(!kinfit_results) {
            const size_t pairId = ntuple::CombinationPairToIndex(selected_bjet_pair, GetNJets());
            const auto iter = std::find(event->kinFit_jetPairId.begin(), event->kinFit_jetPairId.end(), pairId);
            if(iter == event->kinFit_jetPairId.end())
                throw exception("Kinfit information for jet pair (%1%, %2%) is not stored for event %3%.")
                    % selected_bjet_pair.first % selected_bjet_pair.second % eventIdentifier;
            const size_t index = static_cast<size_t>(std::distance(event->kinFit_jetPairId.begin(), iter));
            kinfit_results = std::shared_ptr<kin_fit::FitResults>(new kin_fit::FitResults());
            kinfit_results->convergence = event->kinFit_convergence.at(index);
            kinfit_results->chi2 = event->kinFit_chi2.at(index);
            kinfit_results->probability = TMath::Prob(kinfit_results->chi2, 2);
            kinfit_results->mass = event->kinFit_m.at(index);
        }
        return *kinfit_results;
    }

    LorentzVector GetResonanceMomentum(bool useSVfit, bool addMET)
    {
        if(useSVfit && addMET)
            throw exception("Can't add MET and with SVfit applied.");
        LorentzVector p4 = GetHiggsTTMomentum(useSVfit) + GetHiggsBB().GetMomentum();
        if(addMET)
            p4 += GetMET().GetMomentum();
        return p4;
    }

    double GetMT2()
    {
        if(!mt2.is_initialized()) {
            mt2 = Calculate_MT2(event->p4_1, event->p4_2, GetHiggsBB().GetFirstDaughter().GetMomentum(),
                                               GetHiggsBB().GetSecondDaughter().GetMomentum(), event->pfMET_p4);
        }
        return *mt2;
    }

    const FatJetCandidate* SelectFatJet(double mass_cut, double deltaR_subjet_cut)
    {
        using FatJet = ntuple::TupleFatJet;
        using SubJet = ntuple::TupleSubJet;
        if(!HasBjetPair()) return nullptr;
        for(const FatJetCandidate& fatJet : GetFatJets()) {
            if(fatJet->m(FatJet::MassType::SoftDrop) < mass_cut) continue;
            if(fatJet->subJets().size() < 2) continue;
            std::vector<SubJet> subJets = fatJet->subJets();
            std::sort(subJets.begin(), subJets.end(), [](const SubJet& j1, const SubJet& j2) -> bool {
                return j1.p4().Pt() > j2.p4().Pt(); });
            std::vector<double> deltaR;
            for(size_t n = 0; n < 2; ++n) {
                for(size_t k = 0; k < 2; ++k) {
                    const auto dR = ROOT::Math::VectorUtil::DeltaR(subJets.at(n).p4(),
                                                                   GetHiggsBB().GetDaughterMomentums().at(k));
                    deltaR.push_back(dR);
                }
            }
            if((deltaR.at(0) < deltaR_subjet_cut && deltaR.at(3) < deltaR_subjet_cut)
                    || (deltaR.at(1) < deltaR_subjet_cut && deltaR.at(2) < deltaR_subjet_cut))
                return &fatJet;
        }
        return nullptr;
    }

    void SetMvaScore(double _mva_score) { mva_score = _mva_score; }
    double GetMvaScore() const { return mva_score; }

protected:
    const Event* event;
    const SummaryInfo* summaryInfo;
    TriggerResults triggerResults;

private:
    EventIdentifier eventIdentifier;
    JetPair selected_bjet_pair;
    bool has_bjet_pair;

    std::list<ntuple::TupleJet> tuple_jets;
    std::shared_ptr<JetCollection> jets;
    std::list<ntuple::TupleFatJet> tuple_fatJets;
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

    EventInfo(const Event& _event, const JetPair& _selected_bjet_pair = JetPair(0, 1),
              const SummaryInfo* _summaryInfo = nullptr) :
        EventInfoBase(_event, _selected_bjet_pair, _summaryInfo)
    {
        if(summaryInfo)
            triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors(channel));
    }

    using EventInfoBase::EventInfoBase;

    const FirstLeg& GetFirstLeg()
    {
        if(!leg1) {
            tuple_leg1 = std::shared_ptr<FirstTupleLeg>(new FirstTupleLeg(*event, 1));
            leg1 = std::shared_ptr<FirstLeg>(new FirstLeg(*tuple_leg1, tuple_leg1->iso()));
        }
        return *leg1;
    }

    const SecondLeg& GetSecondLeg()
    {
        if(!leg2) {
            tuple_leg2 = std::shared_ptr<SecondTupleLeg>(new SecondTupleLeg(*event, 2));
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

private:
    std::shared_ptr<FirstTupleLeg> tuple_leg1;
    std::shared_ptr<FirstLeg> leg1;
    std::shared_ptr<SecondTupleLeg> tuple_leg2;
    std::shared_ptr<SecondLeg> leg2;
    std::shared_ptr<HiggsTTCandidate> higgs_tt, higgs_tt_sv;
};

inline std::shared_ptr<EventInfoBase> MakeEventInfo(
        Channel channel, const EventInfoBase::Event& event,
        const EventInfoBase::JetPair& selected_bjet_pair = EventInfoBase::JetPair(0, 1),
        const SummaryInfo* summaryInfo = nullptr)
{
    if(channel == Channel::ETau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<ElectronCandidate, TauCandidate>(
                event, selected_bjet_pair, summaryInfo));
    if(channel == Channel::MuTau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<MuonCandidate, TauCandidate>(
                event, selected_bjet_pair, summaryInfo));
    if(channel == Channel::TauTau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<TauCandidate, TauCandidate>(
                event, selected_bjet_pair, summaryInfo));
    if(channel == Channel::MuMu)
        return std::shared_ptr<EventInfoBase>(new EventInfo<MuonCandidate, MuonCandidate>(
                event, selected_bjet_pair, summaryInfo));
    throw exception("Unsupported channel %1%.") % channel;
}

} // namespace analysis
