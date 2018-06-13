/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <list>
#include <algorithm>
#include "EventTuple.h"
#include "AnalysisTypes.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "KinFitInterface.h"
#include "Candidate.h"
#include "TupleObjects.h"
#include "TriggerResults.h"
#include "SummaryTuple.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "hh-bbtautau/Analysis/include/MT2.h"
#include "h-tautau/McCorrections/include/JECUncertaintiesWrapper.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Cuts/include/Btag_2016.h"

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

    explicit SummaryInfo(const ProdSummary& _summary, const std::string& _uncertainties_source = "")
        : summary(_summary)
    {
        for(size_t n = 0; n < summary.triggers_channel.size(); ++n) {
            const int channel_id = summary.triggers_channel.at(n);
            const Channel channel = static_cast<Channel>(channel_id);
            if(!triggerDescriptors.count(channel))
                triggerDescriptors[channel] = std::make_shared<TriggerDescriptorCollection>();
            triggerDescriptors[channel]->Add(summary.triggers_pattern.at(n), {});
        }
        if(!_uncertainties_source.empty())
            jecUncertainties = std::make_shared<jec::JECUncertaintiesWrapper>(_uncertainties_source);
    }

    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors(Channel channel) const
    {
        if(!triggerDescriptors.count(channel))
            throw exception("Information for channel %1% not found.") % channel;
        return triggerDescriptors.at(channel);
    }

    const ProdSummary& operator*() const { return summary; }
    const ProdSummary* operator->() const { return &summary; }

    const jec::JECUncertaintiesWrapper& GetJecUncertainties() const
    {
        if(!jecUncertainties)
            throw exception("Jec Uncertainties not stored.");
        return *jecUncertainties;
    }

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

        SelectedSignalJets() : selectedBjetPair(ntuple::UndefinedJetPair()),
            selectedVBFjetPair(ntuple::UndefinedJetPair()), n_bjets(0) { }

        bool HasBjetPair(size_t njets) const
        {
            return selectedBjetPair.first < njets && selectedBjetPair.second < njets;
        }

        bool HasVBFPair(size_t njets) const
        {
            return selectedVBFjetPair.first < njets && selectedVBFjetPair.second < njets;
        }

        bool isSelectedBjet(size_t n) const
        {
            return selectedBjetPair.first == n || selectedBjetPair.second == n;
        }

        bool isSelectedVBFjet(size_t n) const
        {
            return selectedVBFjetPair.first == n || selectedVBFjetPair.second == n;
        }

    };

    static SelectedSignalJets SelectSignalJets(const Event& event,
                                               const analysis::Period& period,
                                               JetOrdering jet_ordering = JetOrdering::CSV)
    {
        double bjet_pt_cut = period == analysis::Period::Run2017 ? cuts::btag_2017::pt : cuts::btag_2016::pt;
        double bjet_eta_cut = period == analysis::Period::Run2017 ? cuts::btag_2017::eta : cuts::btag_2016::eta;
        double vbf_pt_cut = 30;
        double vbf_eta_cut = 5;
        double btag_cut;

        if(period == analysis::Period::Run2017){
            btag_cut = jet_ordering == JetOrdering::DeepCSV ? cuts::btag_2017::deepCSVv2M : cuts::btag_2017::CSVv2M ;
        }
        else if(period == analysis::Period::Run2016){
            btag_cut = jet_ordering == JetOrdering::DeepCSV ? cuts::btag_2016::DeepCSVM : cuts::btag_2016::CSVv2M ;
        }
        else
            throw exception("Unsupported btag cut for a specific period.");

        SelectedSignalJets selected_signal_jets;

        const auto CreateJetInfo = [&]() -> auto {
            std::vector<analysis::jet_ordering::JetInfo<decltype(event.jets_p4)::value_type>> jet_info_vector;
            for(size_t n = 0; n < event.jets_p4.size(); ++n) {
                double tag;
                if(jet_ordering == JetOrdering::Pt)
                    tag = event.jets_p4.at(n).Pt();
                else if(jet_ordering == JetOrdering::CSV)
                    tag = event.jets_csv.at(n);
                else if(jet_ordering == JetOrdering::DeepCSV)
                    tag = event.jets_deepCsv_BvsAll.at(n);
                else
                    throw exception("Unsupported jet ordering for b-jet pair selection.");
                jet_info_vector.emplace_back(event.jets_p4.at(n),n,tag);
            }
            return jet_info_vector;
        };


        auto jet_info_vector = CreateJetInfo();
        auto bjets_ordered = jet_ordering::OrderJets(jet_info_vector,true,bjet_pt_cut,bjet_eta_cut);
        selected_signal_jets.n_bjets = bjets_ordered.size();
        if(bjets_ordered.size() >= 1){
            selected_signal_jets.selectedBjetPair.first = bjets_ordered.at(0).index;
        }
        
        if(bjets_ordered.size() >= 2){
            double tag_second_bjet = jet_ordering == JetOrdering::DeepCSV ?
            event.jets_deepCsv_BvsAll.at(bjets_ordered.at(1).index) : event.jets_csv.at(bjets_ordered.at(1).index);
            if(tag_second_bjet > btag_cut){
                selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
            }
        }

        std::vector<analysis::jet_ordering::JetInfo<decltype(event.jets_p4)::value_type>> jet_info_vector_vbf;
        for(size_t n = 0; n < event.jets_p4.size(); ++n) {
            if(selected_signal_jets.isSelectedBjet(n)) continue;
            double tag = event.jets_p4.at(n).Pt();
            jet_info_vector_vbf.emplace_back(event.jets_p4.at(n),n,tag);
        }

        auto vbf_jets_ordered = jet_ordering::OrderJets(jet_info_vector_vbf,true,vbf_pt_cut,vbf_eta_cut);
        
        double max_mjj = -std::numeric_limits<double>::infinity();
        for(size_t n = 0; n < vbf_jets_ordered.size(); ++n) {
            const auto& jet_1 = vbf_jets_ordered.at(n);
            for(size_t h = n+1; h < vbf_jets_ordered.size(); ++h) {
                const auto& jet_2 = vbf_jets_ordered.at(h);
                const auto jet_12 = jet_1.p4 + jet_2.p4;
                if(jet_12.M() > max_mjj){
                    max_mjj = jet_12.M();
                    selected_signal_jets.selectedVBFjetPair = std::make_pair(vbf_jets_ordered.at(n).index,
                                                                             vbf_jets_ordered.at(h).index);
                }
            }
        }
        

        if(selected_signal_jets.HasBjetPair(event.jets_p4.size()) ||
                (!selected_signal_jets.HasBjetPair(event.jets_p4.size()) &&
                 !selected_signal_jets.HasVBFPair(event.jets_p4.size())) ) return selected_signal_jets;

        auto jet_info_vector_new = CreateJetInfo();
        auto new_bjets_ordered = jet_ordering::OrderJets(jet_info_vector_new,true,bjet_pt_cut,bjet_eta_cut);
        if(new_bjets_ordered.size() >= 1)
            selected_signal_jets.selectedBjetPair.second = new_bjets_ordered.at(0).index;
        else{
            selected_signal_jets.selectedVBFjetPair = ntuple::UndefinedJetPair();
            if (bjets_ordered.size() >= 2)
                selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
        }

        return selected_signal_jets;
    }

    std::array<size_t,2> GetSelectedBjetIndices() const
    {
        std::array<size_t,2> bjet_indexes;
        bjet_indexes[0] = selected_signal_jets.selectedBjetPair.first;
        bjet_indexes[1] = selected_signal_jets.selectedBjetPair.second;
        return bjet_indexes;
    }

    std::set<size_t> GetSelectedBjetIndicesSet() const
    {
        std::set<size_t> bjet_indexes;
        bjet_indexes.insert(selected_signal_jets.selectedBjetPair.first);
        bjet_indexes.insert(selected_signal_jets.selectedBjetPair.second);
        return bjet_indexes;
    }


    static constexpr int verbosity = 0;

    EventInfoBase(const Event& _event, Period _period, JetOrdering _jet_ordering,
                  const SummaryInfo* _summaryInfo = nullptr) :
        event(&_event), summaryInfo(_summaryInfo), eventIdentifier(_event.run, _event.lumi, _event.evt),
        selected_signal_jets(SelectSignalJets(_event,_period,_jet_ordering)), period(_period),
        jet_ordering(_jet_ordering)
    {
        mutex = std::make_shared<Mutex>();
        triggerResults.SetAcceptBits(event->trigger_accepts);
        triggerResults.SetMatchBits(event->trigger_matches);
    }
    
    EventInfoBase(const EventInfoBase& ) = default; //copy constructor
    virtual ~EventInfoBase(){} //destructor
    
    EventInfoBase& operator= ( const EventInfoBase& ) = default; //assignment


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

    static const kin_fit::FitProducer& GetKinFitProducer()
    {
        static kin_fit::FitProducer kinfitProducer;
        return kinfitProducer;
    }

    virtual const AnalysisObject& GetLeg(size_t /*leg_id*/) { throw exception("Method not supported."); }
    virtual LorentzVector GetHiggsTTMomentum(bool /*useSVfit*/) { throw exception("Method not supported."); }

    size_t GetNJets() const { return event->jets_p4.size(); }
    size_t GetNFatJets() const { return event->fatJets_p4.size(); }

    const SelectedSignalJets& GetSelectedSignalJets() const {return selected_signal_jets; }
    Period GetPeriod() const { return period; }
    JetOrdering GetJetOrdering() const {return jet_ordering; }

    const JetCollection& GetJets()
    {
        Lock lock(*mutex);
        if(!jets) {
            jets = std::shared_ptr<JetCollection>(new JetCollection());
            tuple_jets = std::make_shared<std::list<ntuple::TupleJet>>();
            for(size_t n = 0; n < GetNJets(); ++n) {
                tuple_jets->push_back(ntuple::TupleJet(*event, n));
                jets->push_back(JetCandidate(tuple_jets->back()));
            }
        }
        return *jets;
    }

    void SetJets(const JetCollection& new_jets)
    {
        Lock lock(*mutex);
        jets = std::make_shared<JetCollection>(new_jets);
    }

    JetCollection SelectJets(double pt_cut = std::numeric_limits<double>::lowest(),
                             double eta_cut = std::numeric_limits<double>::max(),
                             double csv_cut = std::numeric_limits<double>::lowest(),
                             JetOrdering jet_ordering = JetOrdering::CSV,
                             const std::set<size_t>& bjet_indexes = {})
    {
        Lock lock(*mutex);
        const JetCollection& all_jets = GetJets();
        JetCollection selected_jets;
        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        for(size_t n = 0; n < all_jets.size(); ++n) {
            const JetCandidate& jet = all_jets.at(n);
            if(bjet_indexes.count(n) || jet->csv() <= csv_cut) continue;
            double tag;
            if(jet_ordering == JetOrdering::Pt)
                tag = jet.GetMomentum().Pt();
            else if(jet_ordering == JetOrdering::CSV)
                tag = jet->csv();
            else if(jet_ordering == JetOrdering::DeepCSV)
                tag = jet->deepcsv();
            else
                throw exception("Unsupported jet ordering for jet selection.");
            jet_info_vector.emplace_back(jet.GetMomentum(),n,tag);
        }

        auto jets_ordered = jet_ordering::OrderJets(jet_info_vector,true,pt_cut,eta_cut);
        for(size_t h = 0; h < jets_ordered.size(); ++h){
            const JetCandidate& jet = all_jets.at(jets_ordered.at(h).index);
            selected_jets.push_back(jet);
        }
        return selected_jets;
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
        Lock lock(*mutex);
        if(!fatJets) {
            fatJets = std::shared_ptr<FatJetCollection>(new FatJetCollection());
            tuple_fatJets = std::make_shared<std::list<ntuple::TupleFatJet>>();
            for(size_t n = 0; n < GetNFatJets(); ++n) {
                tuple_fatJets->push_back(ntuple::TupleFatJet(*event, n));
                fatJets->push_back(FatJetCandidate(tuple_fatJets->back()));
            }
        }
        return *fatJets;
    }

    bool HasBjetPair() const { return selected_signal_jets.HasBjetPair(GetNJets()); }

    bool HasVBFjetPair() const { return selected_signal_jets.HasVBFPair(GetNJets()); }

    const JetCandidate& GetVBFJet(const size_t index)
    {
        if(!HasVBFjetPair() || (index != 1 && index != 2))
            throw exception("VBF jet not found.");
        if(index == 1)
            return GetJets().at(selected_signal_jets.selectedVBFjetPair.first);
        return GetJets().at(selected_signal_jets.selectedVBFjetPair.second);
    }

    const JetCandidate& GetBJet(const size_t index)
    {
        if(!HasBjetPair() || (index != 1 && index != 2) )
            throw exception("B jet not found.");
        if(index == 1)
            return GetJets().at(selected_signal_jets.selectedBjetPair.first);
        return GetJets().at(selected_signal_jets.selectedBjetPair.second);
    }

    const HiggsBBCandidate& GetHiggsBB()
    {
        Lock lock(*mutex);
        if(!HasBjetPair())
            throw exception("Can't create H->bb candidate.");
        if(!higgs_bb) {
            const auto& jets = GetJets();
            higgs_bb = std::make_shared<HiggsBBCandidate>(jets.at(selected_signal_jets.selectedBjetPair.first),
                                                          jets.at(selected_signal_jets.selectedBjetPair.second));
        }
        return *higgs_bb;
    }

    const MET& GetMET()
    {
        Lock lock(*mutex);
        if(!met) {
            tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
            met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
        }
        return *met;
    }

    template<typename LorentzVector>
    void SetMetMomentum(const LorentzVector& new_met_p4)
    {
        Lock lock(*mutex);
        if(!tuple_met)
            tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
        met = std::make_shared<MET>(*tuple_met, tuple_met->cov());
        met->SetMomentum(new_met_p4);
    }

    const kin_fit::FitResults& GetKinFitResults()
    {
        Lock lock(*mutex);
        if(!HasBjetPair())
            throw exception("Can't retrieve KinFit results.");
        if(!kinfit_results) {
            const size_t pairId = ntuple::CombinationPairToIndex(selected_signal_jets.selectedBjetPair, GetNJets());
            const auto iter = std::find(event->kinFit_jetPairId.begin(), event->kinFit_jetPairId.end(), pairId);
            kinfit_results = std::make_shared<kin_fit::FitResults>();
            if(iter == event->kinFit_jetPairId.end()){
                double energy_resolution_1 = GetBJet(1)->resolution()*GetBJet(1).GetMomentum().E();
                double energy_resolution_2 = GetBJet(2)->resolution()*GetBJet(2).GetMomentum().E();
                const auto& kinfitProducer = GetKinFitProducer();
                const auto& result = kinfitProducer.Fit(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                                                         GetBJet(1).GetMomentum(),
                                                         GetBJet(2).GetMomentum(),
                                                         *met,energy_resolution_1,energy_resolution_2);
                kinfit_results->convergence = result.convergence;
                kinfit_results->chi2 = result.chi2;
                kinfit_results->probability = TMath::Prob(result.chi2, 2);
                kinfit_results->mass = result.mass;
            }
            else {
                const size_t index = static_cast<size_t>(std::distance(event->kinFit_jetPairId.begin(), iter));
                kinfit_results->convergence = event->kinFit_convergence.at(index);
                kinfit_results->chi2 = event->kinFit_chi2.at(index);
                kinfit_results->probability = TMath::Prob(kinfit_results->chi2, 2);
                kinfit_results->mass = event->kinFit_m.at(index);
            }
        }
        return *kinfit_results;
    }

    LorentzVector GetResonanceMomentum(bool useSVfit, bool addMET)
    {
        Lock lock(*mutex);
        if(useSVfit && addMET)
            throw exception("Can't add MET and with SVfit applied.");
        LorentzVector p4 = GetHiggsTTMomentum(useSVfit) + GetHiggsBB().GetMomentum();
        if(addMET)
            p4 += GetMET().GetMomentum();
        return p4;
    }

    double GetMT2()
    {
        Lock lock(*mutex);
        if(!mt2.is_initialized()) {
            mt2 = Calculate_MT2(event->p4_1, event->p4_2, GetHiggsBB().GetFirstDaughter().GetMomentum(),
                                               GetHiggsBB().GetSecondDaughter().GetMomentum(), event->pfMET_p4);
        }
        return *mt2;
    }

    const FatJetCandidate* SelectFatJet(double mass_cut, double deltaR_subjet_cut)
    {
        Lock lock(*mutex);
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

    void SetMvaScore(double _mva_score)
    {
        Lock lock(*mutex);
        mva_score = _mva_score;
    }

    double GetMvaScore() const { return mva_score; }

    virtual std::shared_ptr<EventInfoBase> ApplyShiftBase(UncertaintySource uncertainty_source,
        UncertaintyScale scale) = 0;

protected:
    const Event* event;
    const SummaryInfo* summaryInfo;
    TriggerResults triggerResults;
    std::shared_ptr<Mutex> mutex;

private:
    EventIdentifier eventIdentifier;
    SelectedSignalJets selected_signal_jets;
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
              const SummaryInfo* _summaryInfo = nullptr) :
        EventInfoBase(_event, _period, _jet_ordering, _summaryInfo)
    {
        if(summaryInfo)
            triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors(channel));
    }

    using EventInfoBase::EventInfoBase;

    const FirstLeg& GetFirstLeg()
    {
        Lock lock(*mutex);
        if(!leg1) {
            tuple_leg1 = std::shared_ptr<FirstTupleLeg>(new FirstTupleLeg(*event, 1));
            leg1 = std::shared_ptr<FirstLeg>(new FirstLeg(*tuple_leg1, tuple_leg1->iso()));
        }
        return *leg1;
    }

    const SecondLeg& GetSecondLeg()
    {
        Lock lock(*mutex);
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

inline std::shared_ptr<EventInfoBase> MakeEventInfo(
        Channel channel, const EventInfoBase::Event& event,
        const Period period, JetOrdering jet_ordering,
        const SummaryInfo* summaryInfo = nullptr)
{
    if(channel == Channel::ETau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<ElectronCandidate, TauCandidate>(
                event, period, jet_ordering, summaryInfo));
    if(channel == Channel::MuTau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<MuonCandidate, TauCandidate>(
                event, period, jet_ordering, summaryInfo));
    if(channel == Channel::TauTau)
        return std::shared_ptr<EventInfoBase>(new EventInfo<TauCandidate, TauCandidate>(
                event, period, jet_ordering, summaryInfo));
    if(channel == Channel::MuMu)
        return std::shared_ptr<EventInfoBase>(new EventInfo<MuonCandidate, MuonCandidate>(
                event, period, jet_ordering, summaryInfo));
    throw exception("Unsupported channel %1%.") % channel;
}

} // namespace analysis
