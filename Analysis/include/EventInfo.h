/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <list>
#include "EventTuple.h"
#include "AnalysisTypes.h"
#include "RootExt.h"
#include "KinFitInterface.h"
#include "Candidate.h"
#include "TupleObjects.h"

namespace analysis {

using ElectronCandidate = LeptonCandidate<ntuple::TupleElectron>;
using MuonCandidate = LeptonCandidate<ntuple::TupleMuon>;
using TauCandidate = LeptonCandidate<ntuple::TupleTau>;
using JetCandidate = Candidate<ntuple::TupleJet>;
using MET = MissingET<ntuple::TupleMet>;

class EventInfoBase {
public:
    using Event = ntuple::Event;
    using BjetPair = std::pair<size_t, size_t>;
    using JetCollection = std::vector<JetCandidate>;
    using JetCollectionPtr = std::shared_ptr<JetCollection>;
    using HiggsBBCandidate = CompositCandidate<JetCandidate, JetCandidate>;

    static size_t NumberOfCombinationPairs(size_t n_bjets)
    {
        return n_bjets * (n_bjets - 1) / 2;
    }

    static size_t CombinationPairToIndex(const BjetPair& pair, size_t n_bjets)
    {
        const size_t min = std::min(pair.first, pair.second);
        const size_t max = std::max(pair.first, pair.second);
        if(n_bjets < 2 || min == max || min >= n_bjets || max >= n_bjets)
            throw exception("bad combination pair (%1%, %2%) for n b-jets = %3%.")
                % pair.first % pair.second % n_bjets;
        return max - 1 + min * (2 * n_bjets - 3 - min) / 2;
    }

    static BjetPair CombinationIndexToPair(size_t index, size_t n_bjets)
    {
        if(n_bjets < 2 || index >= n_bjets * (n_bjets - 1) / 2)
            throw exception("bad combination index = %1% for n b-jets = %2%.") % index % n_bjets;

        for(size_t min = 0;; ++min) {
            const size_t l = CombinationPairToIndex(BjetPair(min, n_bjets - 1), n_bjets);
            if(l >= index) {
                const size_t max = index + n_bjets - 1 - l;
                return BjetPair(min, max);
            }
        }
    }

    static constexpr int verbosity = 0;

    EventInfoBase(const Event& _event, const BjetPair& _selected_bjet_pair)
        : event(&_event), selected_bjet_pair(_selected_bjet_pair),
          has_bjet_pair(selected_bjet_pair.first < _event.pt_jets.size() &&
                        selected_bjet_pair.second < _event.pt_jets.size()) {}
    virtual ~EventInfoBase() {}

    const Event& operator*() const { return *event; }
    const Event* operator->() const { return event; }

    EventEnergyScale GetEnergyScale() const { return static_cast<EventEnergyScale>(event->eventEnergyScale); }

    virtual const AnalysisObject& GetLeg(size_t leg_id) = 0;
    virtual LorentzVector GetHiggsTTMomentum(bool useSVfit) = 0;

    const JetCollection& GetJets()
    {
        if(!jets) {
            jets = std::shared_ptr<JetCollection>(new JetCollection());
            for(size_t n = 0; n < event->pt_jets.size(); ++n) {
                tuple_jets.push_back(ntuple::TupleJet(*event, n));
                jets.push_back(JetCandidate(tuple_jets.back()));
            }
        }
        return *jets;
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
        if(!kinfit_resutls) {
            size_t index = CombinationPairToIndex(selected_bjet_pair, event->pt_jets.size());
            kinfit_results = std::shared_ptr<kin_fit::FitResults>(new kin_fit::FitResults());
            kinfit_resutls->convergence = event->kinFit_convergence.at(index);
            kinfit_results->chi2 = event->kinFit_chi2;
            kinfit_results->probability = event->kinFit_probability;
            kinfit_results->mass = event->kinFit_m;
        }
        return *kinfit_resutls;
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
    BjetPair selected_bjet_pair;
    bool has_bjet_pair;

    std::list<ntuple::TupleJet> tuple_jets;
    std::shared_ptr<JetCollection> jets;
    std::shared_ptr<HiggsBBCandidate> higgs_bb;
    std::shared_ptr<ntuple::TupleMet> tuple_met;
    std::shared_ptr<MET> met;
    std::shared_ptr<kin_fit::FitResults> kinfit_resutls;
};

template<typename _FirstLeg>
class EventInfo : public EventInfoBase {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = TauCandidate;
    using FirstTupleLeg = typename FirstLeg::PATObject;
    using SecondTupleLeg = typename SecondLeg::PATObject;
    using HiggsTTCandidate = CompositCandidate<FirstLeg, SecondLeg>;

    const FirstLeg& GetFirstLeg()
    {
        if(!leg1) {
            tuple_leg1 = std::shared_ptr<FirstTupleLeg>(new FirstTupleLeg(*event, 1));
            leg1 = std::shared_ptr<FirstLeg>(new FirstLeg(*tuple_leg1));
        }
        return *leg1;
    }

    const FirstLeg& GetSecondLeg()
    {
        if(!leg2) {
            tuple_leg2 = std::shared_ptr<SecondTupleLeg>(new SecondTupleLeg(*event, 2));
            leg2 = std::shared_ptr<SecondLeg>(new SecondLeg(*tuple_leg2));
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
                const LorentzVectorM p4_sv(event->pt_sv, event->eta_sv, event->phi_sv, event->m_sv);
                higgs_tt_sv = std::shared_ptr<HiggsTTCandidate>(
                              new HiggsTTCandidate(GetFirstLeg(), GetSecondLeg(), p4_sv));
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
