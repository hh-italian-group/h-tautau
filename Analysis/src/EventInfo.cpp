/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventInfo.h"

namespace analysis {

SummaryInfo::SummaryInfo(const ProdSummary& _summary, const std::string& _uncertainties_source)
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

std::shared_ptr<const TriggerDescriptorCollection> SummaryInfo::GetTriggerDescriptors(Channel channel) const
{
    if(!triggerDescriptors.count(channel))
        throw exception("Information for channel %1% not found.") % channel;
    return triggerDescriptors.at(channel);
}

const SummaryInfo::ProdSummary& SummaryInfo::operator*() const { return summary; }
const SummaryInfo::ProdSummary* SummaryInfo::operator->() const { return &summary; }

const jec::JECUncertaintiesWrapper& SummaryInfo::GetJecUncertainties() const
{
    if(!jecUncertainties)
        throw exception("Jec Uncertainties not stored.");
    return *jecUncertainties;
}

EventInfoBase::SelectedSignalJets::SelectedSignalJets() : selectedBjetPair(ntuple::UndefinedJetPair()),
    selectedVBFjetPair(ntuple::UndefinedJetPair()), n_bjets(0) { }

bool EventInfoBase::SelectedSignalJets::HasBjetPair(size_t njets) const
{
    return selectedBjetPair.first < njets && selectedBjetPair.second < njets;
}

bool EventInfoBase::SelectedSignalJets::HasVBFPair(size_t njets) const
{
    return selectedVBFjetPair.first < njets && selectedVBFjetPair.second < njets;
}

bool EventInfoBase::SelectedSignalJets::isSelectedBjet(size_t n) const
{
    return selectedBjetPair.first == n || selectedBjetPair.second == n;
}

bool EventInfoBase::SelectedSignalJets::isSelectedVBFjet(size_t n) const
{
    return selectedVBFjetPair.first == n || selectedVBFjetPair.second == n;
}

EventInfoBase::SelectedSignalJets EventInfoBase::SelectSignalJets(const Event& event,
                                                                  const analysis::Period& period,
                                                                  JetOrdering jet_ordering)
{
    BTagger bTagger(period, jet_ordering);
    const double bjet_pt_cut = bTagger.PtCut();
    const double bjet_eta_cut = bTagger.EtaCut();

    SelectedSignalJets selected_signal_jets;

    const auto CreateJetInfo = [&](bool useBTag) -> auto {
        std::vector<analysis::jet_ordering::JetInfo<decltype(event.jets_p4)::value_type>> jet_info_vector;
        for(size_t n = 0; n < event.jets_p4.size(); ++n) {
            if(selected_signal_jets.isSelectedBjet(n)) continue;
            if(selected_signal_jets.isSelectedVBFjet(n)) continue;
            if(!PassEcalNoiceVetoJets(event.jets_p4.at(n), period, event.jets_pu_id.at(n))) continue;
            if((event.jets_pu_id.at(n) & (1 << 2)) == 0) continue;
//            if(useBTag && (event.jets_pu_id.at(n) & (1 << 2)) == 0) continue;

            const double tag = useBTag ? bTagger.BTag(event,n) : event.jets_p4.at(n).Pt();
            jet_info_vector.emplace_back(event.jets_p4.at(n),n,tag);
        }
        return jet_info_vector;
    };

    auto jet_info_vector = CreateJetInfo(true);
    auto bjets_ordered = jet_ordering::OrderJets(jet_info_vector,true,bjet_pt_cut,bjet_eta_cut);
    selected_signal_jets.n_bjets = bjets_ordered.size();
    if(bjets_ordered.size() >= 1){
        selected_signal_jets.selectedBjetPair.first = bjets_ordered.at(0).index;
    }

    if(bjets_ordered.size() >= 2){
        if(bTagger.Pass(event,bjets_ordered.at(1).index)){
            selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
        }
    }

    auto jet_info_vector_vbf = CreateJetInfo(false);
    auto vbf_jets_ordered = jet_ordering::OrderJets(jet_info_vector_vbf,true,
                                                    cuts::hh_bbtautau_2017::jetID::vbf_pt_cut,
                                                    cuts::hh_bbtautau_2017::jetID::vbf_eta_cut);

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


    if(selected_signal_jets.HasBjetPair(event.jets_p4.size())) return selected_signal_jets;

    auto jet_info_vector_new = CreateJetInfo(true);
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

std::array<size_t,2> EventInfoBase::GetSelectedBjetIndices() const
{
    std::array<size_t,2> bjet_indexes;
    bjet_indexes[0] = selected_signal_jets.selectedBjetPair.first;
    bjet_indexes[1] = selected_signal_jets.selectedBjetPair.second;
    return bjet_indexes;
}

std::set<size_t> EventInfoBase::GetSelectedBjetIndicesSet() const
{
    std::set<size_t> bjet_indexes;
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.first);
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.second);
    return bjet_indexes;
}

const analysis::LepCandidate& EventInfoBase::GetFirstLeg()
{
    Lock lock(*mutex);
    if(!leg1) {
        tuple_leg1 = std::make_shared<ntuple::TupleLepton>(*event, GetLegIndex(1));
        leg1 = std::make_shared<analysis::LepCandidate>(*tuple_leg1, tuple_leg1->iso());
    }
    return *leg1;
}

const analysis::LepCandidate& EventInfoBase::GetSecondLeg()
{
    Lock lock(*mutex);
    if(!leg2) {
        tuple_leg2 = std::make_shared<ntuple::TupleLepton>(*event, GetLegIndex(2));
        leg2 = std::make_shared<analysis::LepCandidate>(*tuple_leg2, tuple_leg2->iso());
    }
    return *leg2;
}

std::shared_ptr<EventInfoBase> EventInfoBase::ApplyShift(UncertaintySource uncertainty_source,
    UncertaintyScale scale)
{
    std::shared_ptr<EventInfoBase> shifted_event_info = std::make_shared<EventInfoBase>(*this);
    const SummaryInfo& summaryInfo = shifted_event_info->GetSummaryInfo();
    const jec::JECUncertaintiesWrapper& jecUncertainties = summaryInfo.GetJecUncertainties();
    const JetCollection& jets = shifted_event_info->GetJets();
    const auto& other_jets_p4 = event->other_jets_p4;
    auto shifted_met_p4(shifted_event_info->GetMET().GetMomentum());
    const JetCollection& corrected_jets = jecUncertainties.ApplyShift(jets,uncertainty_source,scale,&other_jets_p4,&shifted_met_p4);
    shifted_event_info->SetJets(corrected_jets);
    shifted_event_info->SetMetMomentum(shifted_met_p4);
    return shifted_event_info;
}

EventInfoBase::EventInfoBase(const Event& _event, const SummaryInfo* _summaryInfo,
                             size_t _selected_htt_index, const SelectedSignalJets& _selected_signal_jets,
                             Period _period, JetOrdering _jet_ordering) :
event(&_event), summaryInfo(_summaryInfo), selected_htt_index(_selected_htt_index), eventIdentifier(_event.run, _event.lumi, _event.evt),
 selected_signal_jets(_selected_signal_jets), period(_period), jet_ordering(_jet_ordering)
{
    mutex = std::make_shared<Mutex>();
    triggerResults.SetAcceptBits(event->trigger_accepts);
    triggerResults.SetMatchBits(event->trigger_matches.at(selected_htt_index));
    if(summaryInfo)
        triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors(EventInfoBase::GetChannel()));
}


const EventInfoBase::Event& EventInfoBase::operator*() const { return *event; }
const EventInfoBase::Event* EventInfoBase::operator->() const { return event; }

const EventIdentifier& EventInfoBase::GetEventId() const { return eventIdentifier; }
EventEnergyScale EventInfoBase::GetEnergyScale() const
{
    return static_cast<EventEnergyScale>(event->eventEnergyScale);
}

const TriggerResults& EventInfoBase::GetTriggerResults() const { return triggerResults; }
const SummaryInfo& EventInfoBase::GetSummaryInfo() const
{
    if(!summaryInfo)
        throw exception("SummaryInfo was not provided for this event.");
    return *summaryInfo;
}

const kin_fit::FitProducer& EventInfoBase::GetKinFitProducer()
{
    static kin_fit::FitProducer kinfitProducer;
    return kinfitProducer;
}

const sv_fit_ana::FitProducer& EventInfoBase::GetSVFitProducer()
{
    static sv_fit_ana::FitProducer svfitProducer;
    return svfitProducer;
}

// const AnalysisObject& EventInfoBase::GetLeg(size_t /*leg_id*/) { throw exception("Method not supported."); }
// LorentzVector EventInfoBase::GetHiggsTTMomentum(bool /*useSVfit*/) { throw exception("Method not supported."); }

size_t EventInfoBase::GetNJets() const { return event->jets_p4.size(); }
size_t EventInfoBase::GetNFatJets() const { return event->fatJets_p4.size(); }

const EventInfoBase::SelectedSignalJets& EventInfoBase::GetSelectedSignalJets() const { return selected_signal_jets; }
Period EventInfoBase::GetPeriod() const { return period; }
JetOrdering EventInfoBase::GetJetOrdering() const {return jet_ordering; }

const EventInfoBase::JetCollection& EventInfoBase::GetJets()
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

void EventInfoBase::SetJets(const JetCollection& new_jets)
{
    Lock lock(*mutex);
    jets = std::make_shared<JetCollection>(new_jets);
}

EventInfoBase::JetCollection EventInfoBase::SelectJets(double pt_cut, double eta_cut, bool applyPu,
                                                       bool passBtag, JetOrdering jet_ordering,
                                                       const std::set<size_t>& jet_to_exclude_indexes,
                                                       double low_eta_cut)
{
    Lock lock(*mutex);
    BTagger bTagger(period,jet_ordering);
    const JetCollection& all_jets = GetJets();
    JetCollection selected_jets;
    std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
    for(size_t n = 0; n < all_jets.size(); ++n) {
        const JetCandidate& jet = all_jets.at(n);
        if(!PassEcalNoiceVetoJets(jet.GetMomentum(), period, event->jets_pu_id.at(n) )) continue;
        if(jet_to_exclude_indexes.count(n)) continue;
        if(applyPu && (event->jets_pu_id.at(n) & (1 << 2)) == 0) continue;
        if(std::abs(jet.GetMomentum().eta()) < low_eta_cut) continue;
        if(passBtag && !bTagger.Pass(*event,n,DiscriminatorWP::Medium)) continue;

        jet_info_vector.emplace_back(jet.GetMomentum(),n,bTagger.BTag(*event,n));
    }
    auto jets_ordered = jet_ordering::OrderJets(jet_info_vector,true,pt_cut,eta_cut);
    for(size_t h = 0; h < jets_ordered.size(); ++h){
        const JetCandidate& jet = all_jets.at(jets_ordered.at(h).index);
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

double EventInfoBase::GetHT(bool includeHbbJets, bool apply_eta_cut)
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

const EventInfoBase::FatJetCollection& EventInfoBase::GetFatJets()
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

bool EventInfoBase::HasBjetPair() const { return selected_signal_jets.HasBjetPair(GetNJets()); }
bool EventInfoBase::HasVBFjetPair() const { return selected_signal_jets.HasVBFPair(GetNJets()); }

const JetCandidate& EventInfoBase::GetVBFJet(const size_t index)
{
    if(!HasVBFjetPair() || (index != 1 && index != 2))
        throw exception("VBF jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedVBFjetPair.first);
    return GetJets().at(selected_signal_jets.selectedVBFjetPair.second);
}

const JetCandidate& EventInfoBase::GetBJet(const size_t index)
{
    if(!HasBjetPair() || (index != 1 && index != 2) )
        throw exception("B jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedBjetPair.first);
    return GetJets().at(selected_signal_jets.selectedBjetPair.second);
}

const EventInfoBase::HiggsBBCandidate& EventInfoBase::GetHiggsBB()
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

const MET& EventInfoBase::GetMET()
{
    Lock lock(*mutex);
    if(!met) {
        tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
        met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
    }
    return *met;
}

size_t EventInfoBase::GetLegIndex(size_t leg_id)
{
    if(leg_id == 1) return event->first_daughter_indexes.at(selected_htt_index);
    if(leg_id == 2) return event->second_daughter_indexes.at(selected_htt_index);
    throw exception("Invalid leg id = %1%.") % leg_id;
}

bool EventInfoBase::PassDefaultLegSelection(const ntuple::TupleLepton& lepton, Channel channel)
{
    if(lepton.leg_type() == LegType::e || lepton.leg_type() == LegType::mu) return true;
    double pt_cut = channel == Channel::TauTau ? cuts::H_tautau_2016::TauTau::tauID::pt : cuts::H_tautau_2016::MuTau::tauID::pt;
    if(lepton.leg_type() == LegType::tau && lepton.p4().pt() > pt_cut) return true;
    return false;
}

boost::optional<size_t> EventInfoBase::GetHiggsCandidateIndex(const ntuple::Event& event, TauIdDiscriminator discr, double DeltaRmin)
{
    std::vector<ntuple::TupleLepton> lepton_candidates;
    for(size_t n = 0; n < event.lep_p4.size(); ++n)
        lepton_candidates.emplace_back(event, n);

    std::vector<size_t> higgs_candidates;
    const double minDeltaR2 = std::pow(DeltaRmin, 2);
    for(size_t n = 0; n < event.first_daughter_indexes.size(); ++n){
        const auto& first_leg = lepton_candidates.at(event.first_daughter_indexes.at(n));
        const auto& second_leg = lepton_candidates.at(event.second_daughter_indexes.at(n));
        if(ROOT::Math::VectorUtil::DeltaR2(first_leg.p4(), second_leg.p4()) <= minDeltaR2) continue;
        const Channel channel = static_cast<Channel>(event.channelId);
        if(!PassDefaultLegSelection(first_leg,channel)) continue;
        if(!PassDefaultLegSelection(second_leg,channel)) continue;
        higgs_candidates.push_back(n);
    }

    const auto Comparitor = [&](size_t h1, size_t h2) -> bool
    {
        if(h1 == h2) return false;
        for(size_t leg_id = 0; leg_id < 2; ++leg_id) {
            const size_t h1_leg_id = leg_id == 0 ? event.first_daughter_indexes.at(h1) : event.second_daughter_indexes.at(h1);
            const size_t h2_leg_id = leg_id == 0 ? event.first_daughter_indexes.at(h2) : event.second_daughter_indexes.at(h2);
            if(h1_leg_id != h2_leg_id) {
                const auto& h1_leg = lepton_candidates.at(h1_leg_id);
                const auto& h2_leg = lepton_candidates.at(h2_leg_id);
                const int iso_cmp = h1_leg.CompareIsolations(h2_leg, discr);
                if(iso_cmp != 0) return iso_cmp == 1;
                if(h1_leg.p4().pt() != h2_leg.p4().pt())
                    return h1_leg.p4().pt() > h2_leg.p4().pt();
            }
        }

        throw analysis::exception("not found a good criteria for best tau pair");
    };

    if(!higgs_candidates.empty()) return *std::min_element(higgs_candidates.begin(), higgs_candidates.end(), Comparitor);
    return boost::optional<size_t>();
}

const kin_fit::FitResults& EventInfoBase::GetKinFitResults()
{
    Lock lock(*mutex);
    if(!HasBjetPair())
        throw exception("Can't retrieve KinFit results.");
    if(!kinfit_results) {
        const size_t pairId = ntuple::CombinationPairToIndex(selected_signal_jets.selectedBjetPair, GetNJets());
        const auto iter = std::find(event->kinFit_jetPairId.begin(), event->kinFit_jetPairId.end(), pairId);
        kinfit_results = std::make_shared<kin_fit::FitResults>();
        if(iter == event->kinFit_jetPairId.end()){
            double energy_resolution_1 = GetBJet(1)->resolution() * GetBJet(1).GetMomentum().E();
            double energy_resolution_2 = GetBJet(2)->resolution() * GetBJet(2).GetMomentum().E();
            const auto& kinfitProducer = GetKinFitProducer();
            const auto& result = kinfitProducer.Fit(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                                                    GetBJet(1).GetMomentum(), GetBJet(2).GetMomentum(),
                                                    GetMET(), energy_resolution_1, energy_resolution_2);
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

const sv_fit_ana::FitResults& EventInfoBase::GetSVFitResults()
{
    Lock lock(*mutex);
    if(!svfit_results){
        const auto iter = std::find(event->SVfit_Higges_indexes.begin(), event->SVfit_Higges_indexes.end(), selected_htt_index);
        svfit_results = std::make_shared<sv_fit_ana::FitResults>();
        if(iter == event->SVfit_Higges_indexes.end()){
            const auto& svfitProducer = GetSVFitProducer();
            const auto& result = svfitProducer.Fit(GetLeg(1),GetLeg(2),GetMET());
            svfit_results->has_valid_momentum = result.has_valid_momentum;
            svfit_results->momentum = result.momentum;
            svfit_results->momentum_error = result.momentum_error;
            svfit_results->transverseMass = result.transverseMass;
            svfit_results->transverseMass_error = result.transverseMass_error;
        }
        else {
            svfit_results->has_valid_momentum = event->SVfit_is_valid.at(selected_htt_index);
            svfit_results->momentum = event->SVfit_p4.at(selected_htt_index);
            svfit_results->momentum_error = event->SVfit_p4_error.at(selected_htt_index);
            svfit_results->transverseMass = event->SVfit_mt.at(selected_htt_index);
            svfit_results->transverseMass_error = event->SVfit_mt_error.at(selected_htt_index);
        }
    }
    return *svfit_results;
}

LorentzVector EventInfoBase::GetResonanceMomentum(bool useSVfit, bool addMET)
{
    Lock lock(*mutex);
    if(useSVfit && addMET)
        throw exception("Can't add MET and with SVfit applied.");
    LorentzVector p4 = GetHiggsTTMomentum(useSVfit) + GetHiggsBB().GetMomentum();
    if(addMET)
        p4 += GetMET().GetMomentum();
    return p4;
}

double EventInfoBase::GetMT2()
{
    Lock lock(*mutex);
    if(!mt2.is_initialized()) {
        mt2 = Calculate_MT2(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                            GetHiggsBB().GetFirstDaughter().GetMomentum(),
                            GetHiggsBB().GetSecondDaughter().GetMomentum(), event->pfMET_p4);
    }
    return *mt2;
}

const FatJetCandidate* EventInfoBase::SelectFatJet(double mass_cut, double deltaR_subjet_cut)
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

void EventInfoBase::SetMvaScore(double _mva_score)
{
    Lock lock(*mutex);
    mva_score = _mva_score;
}

double EventInfoBase::GetMvaScore() const { return mva_score; }

boost::optional<EventInfoBase> CreateEventInfo(const ntuple::Event& event, const SummaryInfo* summaryInfo,
                                               TauIdDiscriminator discr,
                                               Period period,
                                               JetOrdering jet_ordering)
{
    EventInfoBase::SelectedSignalJets selected_signal_jets  = EventInfoBase::SelectSignalJets(event,period,jet_ordering);
    boost::optional<size_t> selected_higgs_index = EventInfoBase::GetHiggsCandidateIndex(event,discr,cuts::H_tautau_2016::DeltaR_betweenSignalObjects);
    if(selected_higgs_index.is_initialized()) return EventInfoBase(event,summaryInfo,*selected_higgs_index,selected_signal_jets,period,jet_ordering);
    return boost::optional<EventInfoBase>();
}



} // namespace analysis
