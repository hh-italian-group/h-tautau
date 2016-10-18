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
#include "AnalysisTools/Core/include/EventIdentifier.h"

namespace analysis {

using ElectronCandidate = LeptonCandidate<ntuple::TupleElectron>;
using MuonCandidate = LeptonCandidate<ntuple::TupleMuon>;
using TauCandidate = LeptonCandidate<ntuple::TupleTau>;
using JetCandidate = Candidate<ntuple::TupleJet>;
using FatJetCandidate = Candidate<ntuple::TupleFatJet>;
using MET = MissingET<ntuple::TupleMet>;

namespace detail {
template<typename FirstLeg>
constexpr Channel IdentifyChannel();

template<>
inline constexpr Channel IdentifyChannel<ElectronCandidate>() { return Channel::ETau; }

template<>
inline constexpr Channel IdentifyChannel<MuonCandidate>() { return Channel::MuTau; }

template<>
inline constexpr Channel IdentifyChannel<TauCandidate>() { return Channel::TauTau; }
}

struct ChannelInfo {
    template<typename FirstLeg>
    static constexpr Channel IdentifyChannel() { return detail::IdentifyChannel<FirstLeg>(); }
};

enum class JetOrdering { NoOrdering, Pt, CSV };

class EventInfoBase {
public:
    using Event = ntuple::Event;
    using BjetPair = std::pair<size_t, size_t>;
    using JetCollection = std::vector<JetCandidate>;
    using FatJetCollection = std::vector<FatJetCandidate>;
    using HiggsBBCandidate = CompositCandidate<JetCandidate, JetCandidate>;

    static size_t NumberOfCombinationPairs(size_t n_bjets)
    {
        return n_bjets * (n_bjets - 1);
    }

    static size_t CombinationPairToIndex(const BjetPair& pair, size_t n_bjets)
    {
        const size_t min = std::min(pair.first, pair.second);
        const size_t max = std::max(pair.first, pair.second);
        if(n_bjets < 2 || min == max || max >= n_bjets)
            throw exception("bad combination pair (%1%, %2%) for n b-jets = %3%.")
                % pair.first % pair.second % n_bjets;
        const size_t corr = pair.first < pair.second ? -1 : 0;
        return pair.first * (n_bjets - 1) + pair.second + corr;
    }

    static BjetPair CombinationIndexToPair(size_t index, size_t n_bjets)
    {
        if(n_bjets < 2 || index >= NumberOfCombinationPairs(n_bjets))
            throw exception("bad combination index = %1% for n b-jets = %2%.") % index % n_bjets;

        BjetPair pair;
        pair.second = index % (n_bjets - 1);
        pair.first = (index - pair.second) / (n_bjets - 1);
        if(pair.first <= pair.second)
            ++pair.second;
        return pair;
    }

    static BjetPair UndefinedBjetPair()
    {
        static BjetPair pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
        return pair;
    }

    static BjetPair SelectBjetPair(const Event& event, double pt_cut = std::numeric_limits<double>::lowest(),
                                   double eta_cut = std::numeric_limits<double>::lowest(),
                                   JetOrdering jet_ordering = JetOrdering::CSV)
    {
        const auto orderer = [&](size_t j1, size_t j2) -> bool {
            if(jet_ordering == JetOrdering::Pt)
                return event.jets_p4.at(j1).Pt() > event.jets_p4.at(j2).Pt();
            if(jet_ordering == JetOrdering::CSV)
                return event.jets_csv.at(j1) > event.jets_csv.at(j2);
            throw exception("Unsupported jet ordering for b-jet pair selection.");
        };

        std::vector<size_t> indexes;
        for(size_t n = 0; n < event.jets_p4.size(); ++n) {
            if(event.jets_p4.at(n).Pt() > pt_cut && std::abs(event.jets_p4.at(n).eta()) < eta_cut)
                indexes.push_back(n);
        }
        std::sort(indexes.begin(), indexes.end(), orderer);
        BjetPair selected_pair = UndefinedBjetPair();
        if(indexes.size() >= 1)
            selected_pair.first = indexes.at(0);
        if(indexes.size() >= 2)
            selected_pair.second = indexes.at(1);
        return selected_pair;
    }

    static constexpr int verbosity = 0;

    explicit EventInfoBase(const Event& _event, const BjetPair& _selected_bjet_pair = UndefinedBjetPair())
        : event(&_event), eventIdentifier(_event.run, _event.lumi, _event.evt), selected_bjet_pair(_selected_bjet_pair),
          has_bjet_pair(selected_bjet_pair.first < GetNJets() &&
                        selected_bjet_pair.second < GetNJets()) {}
    virtual ~EventInfoBase() {}

    const Event& operator*() const { return *event; }
    const Event* operator->() const { return event; }

    const EventIdentifier& GetEventId() const { return eventIdentifier; }
    EventEnergyScale GetEnergyScale() const { return static_cast<EventEnergyScale>(event->eventEnergyScale); }

    virtual const AnalysisObject& GetLeg(size_t leg_id) { throw exception("Method not supported."); }
    virtual LorentzVector GetHiggsTTMomentum(bool useSVfit) { throw exception("Method not supported."); }

    size_t GetNJets() const { return event->jets_p4.size(); }

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
                             double eta_cut = std::numeric_limits<double>::lowest(),
                             double csv_cut = std::numeric_limits<double>::lowest(),
                             JetOrdering jet_ordering = JetOrdering::CSV)
    {
        const auto orderer = [&](const JetCandidate& j1, const JetCandidate& j2) -> bool {
            if(jet_ordering == JetOrdering::Pt)
                return j1.GetMomentum().Pt() > j2.GetMomentum().Pt();
            return j1->csv() > j2->csv();
        };

        const JetCollection& all_jets = GetJets();
        JetCollection selected_jets;
        for(const JetCandidate& jet : all_jets) {
            if(jet.GetMomentum().Pt() > pt_cut && std::abs(jet.GetMomentum().eta()) < eta_cut
                    && jet->csv() > csv_cut)
                selected_jets.push_back(jet);
        }

        if(jet_ordering != JetOrdering::NoOrdering)
            std::sort(selected_jets.begin(), selected_jets.end(), orderer);
        return selected_jets;
    }

    const FatJetCollection& GetFatJets()
    {
        if(!fatJets) {
            fatJets = std::shared_ptr<FatJetCollection>(new FatJetCollection());
            for(size_t n = 0; n < GetNJets(); ++n) {
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
            size_t index = CombinationPairToIndex(selected_bjet_pair, GetNJets());
            kinfit_results = std::shared_ptr<kin_fit::FitResults>(new kin_fit::FitResults());
            kinfit_results->convergence = event->kinFit_convergence.at(index);
            kinfit_results->chi2 = event->kinFit_chi2.at(index);
            kinfit_results->probability = event->kinFit_probability.at(index);
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

protected:
    const Event* event;

private:
    EventIdentifier eventIdentifier;
    BjetPair selected_bjet_pair;
    bool has_bjet_pair;

    std::list<ntuple::TupleJet> tuple_jets;
    std::shared_ptr<JetCollection> jets;
    std::list<ntuple::TupleFatJet> tuple_fatJets;
    std::shared_ptr<FatJetCollection> fatJets;
    std::shared_ptr<HiggsBBCandidate> higgs_bb;
    std::shared_ptr<ntuple::TupleMet> tuple_met;
    std::shared_ptr<MET> met;
    std::shared_ptr<kin_fit::FitResults> kinfit_results;
};

template<typename _FirstLeg>
class EventInfo : public EventInfoBase {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = TauCandidate;
    using FirstTupleLeg = typename FirstLeg::PATObject;
    using SecondTupleLeg = typename SecondLeg::PATObject;
    using HiggsTTCandidate = CompositCandidate<FirstLeg, SecondLeg>;

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

} // namespace analysis
