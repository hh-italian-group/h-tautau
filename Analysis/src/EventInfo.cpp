/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventInfo.h"

namespace analysis {

SummaryInfo::SummaryInfo(const ProdSummary& _summary, const Channel& _channel,
                         const std::string& _trigger_cfg) :
                         summary(_summary)
{
    if(!_trigger_cfg.empty()){
        triggerDescriptors = TriggerDescriptorCollection::Load(_trigger_cfg,_channel);
    }

}

std::shared_ptr<const TriggerDescriptorCollection> SummaryInfo::GetTriggerDescriptors() const
{
    return triggerDescriptors;
}

const SummaryInfo::ProdSummary& SummaryInfo::operator*() const { return summary; }
const SummaryInfo::ProdSummary* SummaryInfo::operator->() const { return &summary; }

const jec::JECUncertaintiesWrapper& SummaryInfo::GetJecUncertainties() const
{
    if(!jecUncertainties)
        throw exception("Jec Uncertainties not stored.");
    return *jecUncertainties;
}

std::array<size_t,2> EventInfo::GetSelectedBjetIndices() const
{
    std::array<size_t,2> bjet_indexes;
    bjet_indexes[0] = selected_signal_jets.selectedBjetPair.first;
    bjet_indexes[1] = selected_signal_jets.selectedBjetPair.second;
    return bjet_indexes;
}

std::set<size_t> EventInfo::GetSelectedBjetIndicesSet() const
{
    std::set<size_t> bjet_indexes;
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.first);
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.second);
    return bjet_indexes;
}

const analysis::LepCandidate& EventInfo::GetFirstLeg()
{
    Lock lock(*mutex);
    return event_candidate.GetLeptons().at(GetLegIndex(1));
}

const analysis::LepCandidate& EventInfo::GetSecondLeg()
{
    Lock lock(*mutex);
    return event_candidate.GetLeptons().at(GetLegIndex(2));
}

EventInfo::EventInfo(EventCandidate&& _event_candidate, const SummaryInfo* _summaryInfo,
                             size_t _selected_htt_index, const SignalObjectSelector::SelectedSignalJets& _selected_signal_jets,
                             Period _period, JetOrdering _jet_ordering) :
event_candidate(_event_candidate), eventCacheProvider(_event_candidate.GetEvent()),summaryInfo(_summaryInfo), selected_htt_index(_selected_htt_index), eventIdentifier(_event_candidate.GetEvent().run, _event_candidate.GetEvent().lumi, _event_candidate.GetEvent().evt),
 selected_signal_jets(_selected_signal_jets), period(_period), jet_ordering(_jet_ordering)
{
    mutex = std::make_shared<Mutex>();
    triggerResults.SetAcceptBits(event_candidate.GetEvent().trigger_accepts);
    triggerResults.SetMatchBits(event_candidate.GetEvent().trigger_matches.at(selected_htt_index));
    if(summaryInfo)
        triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors());
}


const EventInfo::Event& EventInfo::operator*() const { return event_candidate.GetEvent(); }
const EventInfo::Event* EventInfo::operator->() const { return &(event_candidate.GetEvent()); }

const EventIdentifier& EventInfo::GetEventId() const { return eventIdentifier; }

const TriggerResults& EventInfo::GetTriggerResults() const { return triggerResults; }
const SummaryInfo& EventInfo::GetSummaryInfo() const
{
    if(!summaryInfo)
        throw exception("SummaryInfo was not provided for this event.");
    return *summaryInfo;
}

size_t EventInfo::GetNJets() const
{
    return event_candidate.GetEvent().jets_p4.size();
}

size_t EventInfo::GetNFatJets() const { return event_candidate.GetEvent().fatJets_p4.size(); }
size_t EventInfo::GetHttIndex() const { return selected_htt_index; }
const SignalObjectSelector::SelectedSignalJets& EventInfo::GetSelectedSignalJets() const { return selected_signal_jets; }
Period EventInfo::GetPeriod() const { return period; }
JetOrdering EventInfo::GetJetOrdering() const {return jet_ordering; }


JetCollection EventInfo::SelectJets(double pt_cut, double eta_cut, bool applyPu,
                                                       bool passBtag, JetOrdering jet_ordering,
                                                       const std::set<size_t>& jet_to_exclude_indexes,
                                                       double low_eta_cut,
                                                       analysis::UncertaintySource unc_source,
                                                       analysis::UncertaintyScale unc_scale)
{
    Lock lock(*mutex);
    BTagger bTagger(period,jet_ordering);
    bool base_ordering = jet_ordering != analysis::JetOrdering::HHJetTag;
    const JetCollection& all_jets = GetJets();
    const ntuple::Event& event = event_candidate.GetEvent();
    JetCollection selected_jets;
    std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
    for(size_t n = 0; n < all_jets.size(); ++n) {
        const JetCandidate& jet = all_jets.at(n);
        if(ROOT::Math::VectorUtil::DeltaR(GetLeg(1).GetMomentum(), jet.GetMomentum()) <= cuts::hh_bbtautau_Run2::DeltaR_Lep_Jet) continue;
        if(ROOT::Math::VectorUtil::DeltaR(GetLeg(2).GetMomentum(), jet.GetMomentum()) <= cuts::hh_bbtautau_Run2::DeltaR_Lep_Jet) continue;
        analysis::DiscriminatorIdResults jet_pu_id = jet->GetPuId();
        if(!SignalObjectSelector::PassEcalNoiceVetoJets(jet.GetMomentum(), period, jet_pu_id)) continue;
        if(jet_to_exclude_indexes.count(n)) continue;
        if(applyPu && jet.GetMomentum().pt() < cuts::hh_bbtautau_Run2::jetID::max_pt_veto
                && !(jet_pu_id.Passed(analysis::DiscriminatorWP::Loose))) continue;
        if(std::abs(jet.GetMomentum().eta()) < low_eta_cut) continue;
        if(passBtag && !bTagger.Pass(event,n,unc_source,unc_scale,DiscriminatorWP::Medium)) continue;

        jet_info_vector.emplace_back(jet.GetMomentum(),n,bTagger.BTag(event,n,unc_source,unc_scale,base_ordering));
    }
    auto jets_ordered = jet_ordering::OrderJets(jet_info_vector,true,pt_cut,eta_cut);
    for(size_t h = 0; h < jets_ordered.size(); ++h){
        const JetCandidate& jet = all_jets.at(jets_ordered.at(h).index);
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

double EventInfo::GetHT(bool includeHbbJets, bool apply_eta_cut)
{
    static constexpr double other_jets_min_pt = 20;
    static constexpr double other_jets_max_eta = 4.7;
    static const std::set<size_t> empty_set = {};

    const double eta_cut = apply_eta_cut ? other_jets_max_eta : 5;
    const std::set<size_t>& jets_to_exclude = includeHbbJets ? empty_set : GetSelectedBjetIndicesSet();

    double ht = 0;
    const auto& jets = SelectJets(other_jets_min_pt,eta_cut,true,false,JetOrdering::DeepCSV,jets_to_exclude);
    for(size_t n = 0; n < jets.size(); ++n) {
        const auto& jet = jets.at(n);
        ht += jet.GetMomentum().pt();
    }
    return ht;
}

const FatJetCollection& EventInfo::GetFatJets()
{
    return GetEventCandidate().GetFatJets();
}

bool EventInfo::HasBjetPair() const { return selected_signal_jets.HasBjetPair(GetNJets()); }
bool EventInfo::HasVBFjetPair() const { return selected_signal_jets.HasVBFPair(GetNJets()); }

const JetCandidate& EventInfo::GetVBFJet(const size_t index)
{
    if(!HasVBFjetPair() || (index != 1 && index != 2))
        throw exception("VBF jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedVBFjetPair.first);
    return GetJets().at(selected_signal_jets.selectedVBFjetPair.second);
}

const JetCandidate& EventInfo::GetBJet(const size_t index)
{
    if(!HasBjetPair() || (index != 1 && index != 2) )
        throw exception("B jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedBjetPair.first);
    return GetJets().at(selected_signal_jets.selectedBjetPair.second);
}

const EventInfo::HiggsBBCandidate& EventInfo::GetHiggsBB()
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

size_t EventInfo::GetLegIndex(size_t leg_id)
{
    if(leg_id == 1) return event_candidate.GetEvent().first_daughter_indexes.at(selected_htt_index);
    if(leg_id == 2) return event_candidate.GetEvent().second_daughter_indexes.at(selected_htt_index);
    throw exception("Invalid leg id = %1%.") % leg_id;
}

const kin_fit::FitResults& EventInfo::GetKinFitResults(bool allow_calc, int verbosity)
{
    Lock lock(*mutex);
    if(!HasBjetPair())
        throw exception("Can't retrieve KinFit results.");
    if(!kinfit_results) {
        kinfit_results = std::make_shared<kin_fit::FitResults>();
        LegPair selected_htt_pair = ntuple::LegIndexToPair(selected_htt_index);
        bool gotKinFit = eventCacheProvider.TryGetKinFit(*kinfit_results,selected_htt_pair,
                                        selected_signal_jets.selectedBjetPair,
                                        event_candidate.GetUncSource(),event_candidate.GetScale());
        if(!allow_calc && !gotKinFit)
            throw exception("Not allowed to calculate KinFit.");
        else if(!gotKinFit){
            double energy_resolution_1 = GetBJet(1)->resolution() * GetBJet(1).GetMomentum().E();
            double energy_resolution_2 = GetBJet(2)->resolution() * GetBJet(2).GetMomentum().E();
            const auto& result = kin_fit::FitProducer::Fit(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                                                           GetBJet(1).GetMomentum(), GetBJet(2).GetMomentum(),
                                                           event_candidate.GetMET(), energy_resolution_1,
                                                           energy_resolution_2, verbosity);
            kinfit_results->convergence = result.convergence;
            kinfit_results->chi2 = result.chi2;
            kinfit_results->probability = TMath::Prob(result.chi2, 2);
            kinfit_results->mass = result.mass;
        }
        else{
            kinfit_results->probability = TMath::Prob(kinfit_results->chi2, 2);
        }
    }
    return *kinfit_results;
}

const sv_fit_ana::FitResults& EventInfo::GetSVFitResults(bool allow_calc, int verbosity)
{
    Lock lock(*mutex);
    if(!svfit_results){
        svfit_results = std::make_shared<sv_fit_ana::FitResults>();
        LegPair selected_htt_pair = ntuple::LegIndexToPair(selected_htt_index);
        bool gotSVFit = eventCacheProvider.TryGetSVFit(*svfit_results,selected_htt_pair,
                                        event_candidate.GetUncSource(),event_candidate.GetScale());
        if(!allow_calc && !gotSVFit)
            throw exception("Not allowed to calculate SVFit.");
        else if(!gotSVFit){
            const auto& result = sv_fit_ana::FitProducer::Fit(GetLeg(1), GetLeg(2), event_candidate.GetMET(),
                                                              verbosity);
            svfit_results->has_valid_momentum = result.has_valid_momentum;
            svfit_results->momentum = result.momentum;
            svfit_results->momentum_error = result.momentum_error;
            svfit_results->transverseMass = result.transverseMass;
            svfit_results->transverseMass_error = result.transverseMass_error;
        }
    }
    return *svfit_results;
}

LorentzVector EventInfo::GetResonanceMomentum(bool useSVfit, bool addMET, bool allow_calc)
{
    Lock lock(*mutex);
    if(useSVfit && addMET)
        throw exception("Can't add MET and with SVfit applied.");
    LorentzVector p4 (0,0,0,0);
    if(!useSVfit || (useSVfit && GetSVFitResults(allow_calc).has_valid_momentum))
        p4 = GetHiggsTTMomentum(useSVfit, allow_calc) + GetHiggsBB().GetMomentum();
    if(addMET)
        p4 += event_candidate.GetMET().GetMomentum();
    return p4;
}

double EventInfo::GetMT2()
{
    Lock lock(*mutex);
    if(!mt2.is_initialized()) {
        mt2 = Calculate_MT2(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                            GetHiggsBB().GetFirstDaughter().GetMomentum(),
                            GetHiggsBB().GetSecondDaughter().GetMomentum(), GetMET().GetMomentum());
    }
    return *mt2;
}

const FatJetCandidate* EventInfo::SelectFatJet(double mass_cut, double deltaR_subjet_cut)
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

void EventInfo::SetMvaScore(double _mva_score)
{
    Lock lock(*mutex);
    mva_score = _mva_score;
}

double EventInfo::GetMvaScore() const { return mva_score; }

const JetCollection& EventInfo::GetJets() { return GetEventCandidate().GetJets(); }
const MET& EventInfo::GetMET() { return GetEventCandidate().GetMET(); }

boost::optional<EventInfo> CreateEventInfo(const ntuple::Event& event,
                                               const SignalObjectSelector& signalObjectSelector,
                                               const SummaryInfo* summaryInfo,
                                               Period period,
                                               JetOrdering jet_ordering,
                                               bool is_sync,
                                               UncertaintySource uncertainty_source,
                                               UncertaintyScale scale,
                                               bool debug)
{
    if(debug)
        std::cout << "Creating event candidate... ";
    EventCandidate event_candidate(event, uncertainty_source, scale);
    if(debug)
        std::cout << "done\nSelecting higgs->tautau candidate...\n";
    boost::optional<size_t> selected_higgs_index =
            signalObjectSelector.GetHiggsCandidateIndex(event_candidate, is_sync);
    if(debug) {
        if(selected_higgs_index.is_initialized())
            std::cout << "Higgs->tautau candidate number " << *selected_higgs_index << " is selected.\n";
        else
            std::cout << "Higgs->tautau candidate is not selected.\n";
    }
    if(!selected_higgs_index.is_initialized()) return boost::optional<EventInfo>();
    if(debug)
        std::cout << "Selecting signal jets... ";
    auto selected_signal_jets  = signalObjectSelector.SelectSignalJets(event_candidate, period, jet_ordering,
                                                                       *selected_higgs_index, uncertainty_source,
                                                                       scale);
    if(debug)
        std::cout << "done\nCreating EventInfo...";
    auto eventInfo = EventInfo(std::move(event_candidate), summaryInfo, *selected_higgs_index,
                                   selected_signal_jets, period, jet_ordering);
    if(debug)
        std::cout << "done\n";
    return eventInfo;
}

} // namespace analysis
