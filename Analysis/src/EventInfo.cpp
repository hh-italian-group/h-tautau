/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Cuts/include/btag_Run2.h"

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

namespace {
    TriggerResults _InitializeTriggerResults(const ntuple::Event& event,
                                             const std::shared_ptr<const SummaryInfo>& summaryInfo,
                                             size_t selected_htt_index)
    {
        TriggerResults triggerResults;
        triggerResults.SetAcceptBits(event.trigger_accepts);
        triggerResults.SetMatchBits(event.trigger_matches.at(selected_htt_index));
        if(summaryInfo)
            triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors());
        return triggerResults;
    }

    ntuple::LegPair _InitializeSelectedTauIndices(const ntuple::Event& event, size_t selected_htt_index)
    {
        if(selected_htt_index >= event.first_daughter_indexes.size()
                || selected_htt_index >= event.second_daughter_indexes.size() )
            throw exception("Tau tau indices for htt_index=%1% not found in ntuple::Event.") % selected_htt_index;
        return ntuple::LegPair(event.first_daughter_indexes.at(selected_htt_index),
                               event.second_daughter_indexes.at(selected_htt_index));
    }
}

EventInfo::EventInfo(const std::shared_ptr<EventCandidate>& _event_candidate,
                     const std::shared_ptr<const SummaryInfo>& _summaryInfo, size_t _selected_htt_index,
                     const SelectedSignalJets& _selected_signal_jets, const BTagger& _bTagger) :
        event_candidate(_event_candidate), summaryInfo(_summaryInfo),
        triggerResults(_InitializeTriggerResults(_event_candidate->GetEvent(), _summaryInfo, _selected_htt_index)),
        bTagger(_bTagger), selected_htt_index(_selected_htt_index),
        selected_tau_indices(_InitializeSelectedTauIndices(_event_candidate->GetEvent(), _selected_htt_index)),
        selected_signal_jets(_selected_signal_jets)
{
}

const EventInfo::Event& EventInfo::operator*() const { return event_candidate->GetEvent(); }
const EventInfo::Event* EventInfo::operator->() const { return &event_candidate->GetEvent(); }

const EventCandidate& EventInfo::GetEventCandidate() const { return *event_candidate; }
const SummaryInfo& EventInfo::GetSummaryInfo() const
{
    if(!summaryInfo)
        throw exception("SummaryInfo was not provided for this event.");
    return *summaryInfo;
}
size_t EventInfo::GetHttIndex() const { return selected_htt_index; }
const ntuple::LegPair& EventInfo::GetSelectedTauIndices() const { return selected_tau_indices; }
const EventInfo::SelectedSignalJets& EventInfo::GetSelectedSignalJets() const { return selected_signal_jets; }
const BTagger& EventInfo::GetBTagger() const { return bTagger; }
const TriggerResults& EventInfo::GetTriggerResults() const { return triggerResults; }

Channel EventInfo::GetChannel() const { return event_candidate->GetChannel(); }
Period EventInfo::GetPeriod() const { return event_candidate->GetPeriod(); }
const EventIdentifier& EventInfo::GetEventId() const { return event_candidate->GetEventId(); }
const MET& EventInfo::GetMET() const { return event_candidate->GetMET(); }

const LepCandidate& EventInfo::GetLeg(size_t leg_id) const
{
    const size_t leg_index = selected_tau_indices.Get(leg_id);
    return event_candidate->GetLeptons().at(leg_index);
}
const LepCandidate& EventInfo::GetFirstLeg() const { return GetLeg(1); }
const LepCandidate& EventInfo::GetSecondLeg() const { return GetLeg(2); }
const EventInfo::HiggsTTCandidate& EventInfo::GetHiggsTT(bool useSVfit, bool allow_calc)
{
    Lock lock(mutex);
    if(useSVfit) {
        if(!higgs_tt_sv) {
            if(!GetSVFitResults(allow_calc).has_valid_momentum) throw exception("SVFit not converged");
            higgs_tt_sv = HiggsTTCandidate(GetFirstLeg(), GetSecondLeg(), GetSVFitResults(allow_calc).momentum);
        }
        return *higgs_tt_sv;
    }
    if(!higgs_tt)
        higgs_tt = HiggsTTCandidate(GetFirstLeg(), GetSecondLeg());
    return *higgs_tt;
}
boost::optional<LorentzVector> EventInfo::GetHiggsTTMomentum(bool useSVfit, bool allow_calc)
{
    boost::optional<LorentzVector> p4;
    if(!useSVfit || GetSVFitResults(allow_calc).has_valid_momentum)
        p4 = GetHiggsTT(useSVfit, allow_calc).GetMomentum();
    return p4;
}

bool EventInfo::HasBjetPair() const { return selected_signal_jets.HasBjetPair(); }
const JetCandidate& EventInfo::GetBJet(size_t index) const
{
    if(!HasBjetPair() || (index != 1 && index != 2) )
        throw exception("B jet not found.");
    return event_candidate->GetJets().at(selected_signal_jets.bjet_pair.Get(index));
}
const EventInfo::HiggsBBCandidate& EventInfo::GetHiggsBB()
{
    Lock lock(mutex);
    if(!higgs_bb) {
        if(!HasBjetPair())
            throw exception("Can't create H->bb candidate.");
        higgs_bb = HiggsBBCandidate(GetBJet(1), GetBJet(2));
    }
    return *higgs_bb;
}

bool EventInfo::HasVBFjetPair() const { return selected_signal_jets.HasVBFPair(); }
const JetCandidate& EventInfo::GetVBFJet(size_t index) const
{
    if(!HasVBFjetPair() || (index != 1 && index != 2))
        throw exception("VBF jet not found.");
    return event_candidate->GetJets().at(selected_signal_jets.vbf_pair.Get(index));
}

boost::optional<LorentzVector> EventInfo::GetResonanceMomentum(bool useSVfit, bool addMET, bool allow_calc)
{
    Lock lock(mutex);
    if(useSVfit && addMET)
        throw exception("Can't add MET and with SVfit applied.");
    boost::optional<LorentzVector> p4;
    if(!useSVfit || GetSVFitResults(allow_calc).has_valid_momentum)
        p4 = *GetHiggsTTMomentum(useSVfit, allow_calc) + GetHiggsBB().GetMomentum();
    if(p4 && addMET)
        p4 = *p4 + event_candidate->GetMET().GetMomentum();
    return p4;
}

double EventInfo::GetHT(bool includeHbbJets, bool includeVBFJets) const
{
    static const auto pt_cut = Cut1D_Bound::L(cuts::hh_bbtautau_Run2::jetID::pt);
    static const auto eta_cut = Cut1D_Bound::AbsU(cuts::hh_bbtautau_Run2::jetID::vbf_eta);

    SelectedSignalJets selected_jets = selected_signal_jets;
    if(includeHbbJets)
        selected_jets.bjet_pair = ntuple::LegPair::Undefined;
    if(includeVBFJets)
        selected_jets.vbf_pair = ntuple::LegPair::Undefined;
    const auto jet_infos = SignalObjectSelector::CreateJetInfos(*event_candidate, bTagger, true,
                                                                selected_htt_index, selected_jets);
    const auto filtered_jets = jet_ordering::FilterJets(jet_infos, pt_cut, eta_cut);
    double ht = 0;
    for(const auto& jet : filtered_jets)
        ht += event_candidate->GetJets().at(jet.index).GetMomentum().pt();
    return ht;
}

const sv_fit_ana::FitResults& EventInfo::GetSVFitResults(bool allow_calc, int verbosity)
{
    Lock lock(mutex);
    if(!svfit_results) {
        svfit_results = sv_fit_ana::FitResults();
        svfit_results = event_candidate->GetCacheProvider().TryGetSVFit(selected_htt_index,
                                                                        event_candidate->GetUncSource(),
                                                                        event_candidate->GetUncScale());
        if(!svfit_results) {
            if(!allow_calc)
                throw exception("Not allowed to calculate SVFit.");
            svfit_results = sv_fit_ana::FitProducer::Fit(GetLeg(1), GetLeg(2), event_candidate->GetMET(), verbosity);
        }
    }
    return *svfit_results;
}

const kin_fit::FitResults& EventInfo::GetKinFitResults(bool allow_calc, int verbosity)
{
    Lock lock(mutex);
    if(!kinfit_results) {
        if(!HasBjetPair())
            throw exception("Can't retrieve KinFit results.");

        kinfit_results = event_candidate->GetCacheProvider().TryGetKinFit(selected_htt_index,
                                                                          selected_signal_jets.bjet_pair.ToIndex(),
                                                                          event_candidate->GetUncSource(),
                                                                          event_candidate->GetUncScale());
        if(!kinfit_results) {
            if(!allow_calc)
                throw exception("Not allowed to calculate KinFit.");
            const double energy_resolution_1 = GetBJet(1)->resolution() * GetBJet(1).GetMomentum().E();
            const double energy_resolution_2 = GetBJet(2)->resolution() * GetBJet(2).GetMomentum().E();
            kinfit_results = kin_fit::FitProducer::Fit(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                                                       GetBJet(1).GetMomentum(), GetBJet(2).GetMomentum(),
                                                       event_candidate->GetMET(), energy_resolution_1,
                                                       energy_resolution_2, verbosity);
        }
        kinfit_results->probability = TMath::Prob(kinfit_results->chi2, 2);
    }
    return *kinfit_results;
}

double EventInfo::GetMT2()
{
    Lock lock(mutex);
    if(!mt2) {
        mt2 = Calculate_MT2(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                            GetHiggsBB().GetFirstDaughter().GetMomentum(),
                            GetHiggsBB().GetSecondDaughter().GetMomentum(), GetMET().GetMomentum());
    }
    return *mt2;
}

const std::vector<const JetCandidate*>& EventInfo::GetCentralJets()
{
    static const auto pt_cut = Cut1D_Bound::L(cuts::btag_Run2::pt);
    static const auto eta_cut = Cut1D_Bound::AbsU(cuts::btag_Run2::eta);
    Lock lock(mutex);
    if(!central_jets) {
        const auto jet_infos = SignalObjectSelector::CreateJetInfos(*event_candidate, bTagger, true, selected_htt_index,
                                                                    SignalObjectSelector::SelectedSignalJets());

        const auto ordered_jet_infos = jet_ordering::OrderJets(jet_infos, true, pt_cut, eta_cut);
        central_jets = std::vector<const JetCandidate*>();
        for(const auto& jet_info : ordered_jet_infos)
            central_jets->push_back(&event_candidate->GetJets().at(jet_info.index));
    }
    return *central_jets;
}

const std::vector<const JetCandidate*>& EventInfo::GetForwardJets()
{
    static const auto pt_cut = Cut1D_Bound::L(cuts::btag_Run2::pt);
    static const Cut1D_Interval eta_cut(Cut1D_Bound::AbsL(cuts::btag_Run2::eta, true),
                                        Cut1D_Bound::AbsU(cuts::hh_bbtautau_Run2::jetID::vbf_eta));
    Lock lock(mutex);
    if(!forward_jets) {
        const BTagger ptTagger(event_candidate->GetPeriod(), BTaggerKind::Pt);
        const auto jet_infos = SignalObjectSelector::CreateJetInfos(*event_candidate, ptTagger, true,
                                                                    selected_htt_index,
                                                                    SignalObjectSelector::SelectedSignalJets());

        const auto ordered_jet_infos = jet_ordering::OrderJets(jet_infos, true, pt_cut, eta_cut);
        forward_jets = std::vector<const JetCandidate*>();
        for(const auto& jet_info : ordered_jet_infos)
            forward_jets->push_back(&event_candidate->GetJets().at(jet_info.index));
    }
    return *forward_jets;
}
const std::vector<const JetCandidate*>& EventInfo::GetAllJets()
{
    static const auto pt_cut = Cut1D_Bound::L(cuts::btag_Run2::pt);
    static const auto eta_cut = Cut1D_Bound::AbsU(cuts::hh_bbtautau_Run2::jetID::vbf_eta);

    Lock lock(mutex);
    if(!all_jets) {
        const BTagger ptTagger(event_candidate->GetPeriod(), BTaggerKind::Pt);
        const auto jet_infos = SignalObjectSelector::CreateJetInfos(*event_candidate, ptTagger, true,
                                                                    selected_htt_index,
                                                                    SignalObjectSelector::SelectedSignalJets());

        const auto ordered_jet_infos = jet_ordering::OrderJets(jet_infos, true, pt_cut, eta_cut);
        all_jets = std::vector<const JetCandidate*>();
        for(const auto& jet_info : ordered_jet_infos)
            all_jets->push_back(&event_candidate->GetJets().at(jet_info.index));
    }
    return *all_jets;
}

void EventInfo::SetMvaScore(double _mva_score)
{
    Lock lock(mutex);
    mva_score = _mva_score;
}

double EventInfo::GetMvaScore() const
{
    if(!mva_score)
        throw exception("EventInfo: mva_score is not set.");
    return *mva_score;
}

std::unique_ptr<EventInfo> EventInfo::Create(const ntuple::Event& event,
                                             const SignalObjectSelector& signalObjectSelector,
                                             const BTagger& bTagger, DiscriminatorWP btag_wp,
                                             SummaryInfoPtr summaryInfo, UncertaintySource unc_source,
                                             UncertaintyScale unc_scale, bool is_sync, bool debug)
{
    if(debug)
        std::cout << "Creating event candidate... ";
    auto event_candidate = std::make_shared<EventCandidate>(event, unc_source, unc_scale);
    return Create(event_candidate, signalObjectSelector, bTagger, btag_wp, summaryInfo, is_sync, debug);
}

std::unique_ptr<EventInfo> EventInfo::Create(const std::shared_ptr<EventCandidate>& event_candidate,
                                             const SignalObjectSelector& signalObjectSelector,
                                             const BTagger& bTagger, DiscriminatorWP btag_wp,
                                             SummaryInfoPtr summaryInfo, bool is_sync, bool debug)
{
    if(debug)
        std::cout << "done\nSelecting higgs->tautau candidate...\n";
    auto htt_index = signalObjectSelector.GetHiggsCandidateIndex(*event_candidate, is_sync);
    if(debug) {
        if(htt_index.is_initialized())
            std::cout << "Higgs->tautau candidate number " << *htt_index << " is selected.\n";
        else
            std::cout << "Higgs->tautau candidate is not selected.\n";
    }
    if(!htt_index.is_initialized()) return std::unique_ptr<EventInfo>();
    if(debug)
        std::cout << "Selecting signal jets... ";
    event_candidate->SetHHTagScores(*htt_index);
    auto selected_signal_jets = signalObjectSelector.SelectSignalJets(*event_candidate, *htt_index, bTagger, btag_wp);
    if(debug)
        std::cout << "done\nCreating EventInfo...";
    auto eventInfo = std::make_unique<EventInfo>(event_candidate, summaryInfo, *htt_index, selected_signal_jets,
                                                 bTagger);
    if(debug)
        std::cout << "done\n";
    return eventInfo;
}

} // namespace analysis
