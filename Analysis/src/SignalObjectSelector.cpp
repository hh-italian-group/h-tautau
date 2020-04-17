/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Cuts/include/H_tautau_Run2.h"

namespace analysis {

namespace jet_ordering {

bool CompareJets(const JetInfo& jet_1, const JetInfo& jet_2, const Cut1D& pt_cut, const Cut1D& eta_cut)
{
    const bool pass_eta_1 = eta_cut(jet_1.p4.eta()), pass_eta_2 = eta_cut(jet_2.p4.eta());
    if(pass_eta_1 && !pass_eta_2) return true;
    if(!pass_eta_1 && pass_eta_2) return false;
    const bool pass_pt_1 = pt_cut(jet_1.p4.pt()), pass_pt_2 = pt_cut(jet_2.p4.pt());
    if(pass_pt_1 && !pass_pt_2) return true;
    if(!pass_pt_1 && pass_pt_2) return false;

    if(jet_1.tag != jet_2.tag)
        return jet_1.tag > jet_2.tag;
    return jet_1.p4.pt() > jet_2.p4.pt();
};

std::vector<JetInfo> FilterJets(const std::vector<JetInfo>& jets, const Cut1D& pt_cut, const Cut1D& eta_cut)
{
    std::vector<JetInfo> filtered_jets;
    for(const JetInfo& jet : jets) {
        if(pt_cut(jet.p4.pt()) && eta_cut(jet.p4.eta()))
            filtered_jets.push_back(jet);
    }
    return filtered_jets;
}

std::vector<JetInfo> OrderJets(const std::vector<JetInfo>& jets, bool apply_hard_cut,
                               const Cut1D& pt_cut, const Cut1D& eta_cut)
{
    const auto comparitor = [&](const JetInfo& jet_1, const JetInfo& jet_2) -> bool {
        return CompareJets(jet_1, jet_2, pt_cut, eta_cut);
    };

    std::vector<JetInfo> ordered_jets;
    if(apply_hard_cut)
        ordered_jets = FilterJets(jets, pt_cut, eta_cut);
    else
        ordered_jets = jets;
    std::sort(ordered_jets.begin(), ordered_jets.end(), comparitor);
    return ordered_jets;
}

} // namespace jet_ordering

SignalObjectSelector::SignalObjectSelector(SignalMode _mode) : mode(_mode)
{
    if(mode == SignalMode::HTT || mode == SignalMode::TauPOG)
        DR2_leptons = std::pow(cuts::H_tautau_Run2::DeltaR_Lep_Lep, 2);
    else if(mode == SignalMode::HH || mode == SignalMode::HH_legacy)
        DR2_leptons = std::pow(cuts::hh_bbtautau_Run2::DeltaR_Lep_Lep, 2);
    else
        throw exception("Signal Mode for SignalObjectSelector constructor not supported");
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSjetDiscriminator() const
{
    static const std::map<SignalMode, std::pair<TauIdDiscriminator, DiscriminatorWP>> tauVSjetDiscriminator = {
        { SignalMode::HTT, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium } },
        { SignalMode::HH, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium } },
        { SignalMode::HH_legacy, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,
                                   DiscriminatorWP::Medium } },
    };

    if(!tauVSjetDiscriminator.count(mode))
        throw exception("Mode'%1%' isn't found at the tauVSjetDiscriminator map") % mode;
    return tauVSjetDiscriminator.at(mode);
}

std::pair<DiscriminatorWP, DiscriminatorWP> SignalObjectSelector::GetTauVSjetSidebandWPRange() const
{
    static const std::map<SignalMode, std::pair<DiscriminatorWP, DiscriminatorWP>> tauVSjetWps = {
        { SignalMode::HTT, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG, { DiscriminatorWP::VVVLoose, DiscriminatorWP::Medium } },
        { SignalMode::HH, { DiscriminatorWP::VVVLoose, DiscriminatorWP::Medium } },
        { SignalMode::HH_legacy, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
    };

    if(!tauVSjetWps.count(mode))
        throw exception("Mode'%1%' isn't found at the tauVSjetWps map") % mode;
    return tauVSjetWps.at(mode);
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSeDiscriminator(Channel channel) const
{
    using Map = std::map<std::pair<Channel, SignalMode>, std::pair<TauIdDiscriminator, DiscriminatorWP>>;
    static const Map tauVSeDiscriminators = {
        { {Channel::ETau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::ETau, SignalMode::TauPOG}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::TauPOG}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} },
        { {Channel::TauTau, SignalMode::TauPOG},
            {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} },
        { {Channel::ETau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVLoose} },
        { {Channel::MuMu, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} },
        { {Channel::ETau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HH_legacy},
            {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
    };

    auto disc_wp = tauVSeDiscriminators.find({channel, mode});
    if(disc_wp == tauVSeDiscriminators.end())
        throw exception("Channel: %1% and mode: %2% not found in the tauVSeDiscriminators map") % channel % mode;
    return disc_wp->second;
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSmuDiscriminator(Channel channel) const
{
    using Map = std::map<std::pair<Channel, SignalMode>, std::pair<TauIdDiscriminator, DiscriminatorWP>>;
    static const Map tauVSmuDiscriminators = {
        { {Channel::ETau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::TauPOG}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::TauPOG}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::TauPOG}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::VLoose} },
        { {Channel::MuMu, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::VLoose} },
        { {Channel::ETau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
    };
    auto disc_wp = tauVSmuDiscriminators.find({channel, mode});
    if(disc_wp == tauVSmuDiscriminators.end())
        throw exception("Channel: %1% and mode: %2% not found in the tauVSmuDiscriminators map") % channel % mode;
    return disc_wp->second;
}

bool SignalObjectSelector::PassLeptonSelection(const LepCandidate& lepton, Channel channel, const size_t legId,
                                               bool is_sync) const
{
    if(mode == SignalMode::HTT)
        return PassHTT_LeptonSelection(lepton, channel, is_sync);
    if(mode == SignalMode::TauPOG)
        return PassTauPOG_LeptonSelection(lepton, channel);
    if(mode == SignalMode::HH || mode == SignalMode::HH_legacy)
        return PassHH_LeptonSelection(lepton, channel, legId, is_sync);
    throw exception("Signal Mode for SignalObjectSelector class not supported");
}

boost::optional<size_t> SignalObjectSelector::GetHiggsCandidateIndex(const EventCandidate& event_candidate,
                                                                     bool is_sync) const
{
    const ntuple::Event& event = event_candidate.GetEvent();
    std::vector<size_t> higgs_candidates;
    for(size_t n = 0; n < event.first_daughter_indexes.size(); ++n){
        if(event.first_daughter_indexes.at(n) >=event_candidate.GetLeptons().size() ||
            event.second_daughter_indexes.at(n) >=event_candidate.GetLeptons().size())
            throw exception("Daughter indexes greater than lepton size.");
        auto& first_leg = event_candidate.GetLeptons().at(event.first_daughter_indexes.at(n));
        auto& second_leg = event_candidate.GetLeptons().at(event.second_daughter_indexes.at(n));
        if(ROOT::Math::VectorUtil::DeltaR2(first_leg.GetMomentum(), second_leg.GetMomentum()) <= DR2_leptons) continue;
        const Channel channel = static_cast<Channel>(event.channelId);
        if(!PassLeptonSelection(first_leg,channel,1, is_sync)) continue;
        if(!PassLeptonSelection(second_leg,channel,2, is_sync)) continue;
        higgs_candidates.push_back(n);
    }
    const auto Comparitor = [&](size_t h1, size_t h2) -> bool
    {
        bool are_identical = true;
        if(h1 == h2) return false;
        for(size_t leg_id = 0; leg_id < 2; ++leg_id) {
            const size_t h1_leg_id = leg_id == 0 ? event.first_daughter_indexes.at(h1) : event.second_daughter_indexes.at(h1);
            const size_t h2_leg_id = leg_id == 0 ? event.first_daughter_indexes.at(h2) : event.second_daughter_indexes.at(h2);

            if(h1_leg_id != h2_leg_id) {
                are_identical = false;
                auto& h1_leg = event_candidate.GetLeptons().at(h1_leg_id);
                auto& h2_leg = event_candidate.GetLeptons().at(h2_leg_id);
                const int iso_cmp = h1_leg->CompareIsolations(*h2_leg, GetTauVSjetDiscriminator().first);
                if(iso_cmp != 0) return iso_cmp == 1;
                if(h1_leg.GetMomentum().pt() != h2_leg.GetMomentum().pt())
                    return h1_leg.GetMomentum().pt() > h2_leg.GetMomentum().pt();
            }
        }
        if(are_identical) return false;
        throw exception("not found a good criteria for best tau pair for %1%") % EventIdentifier(event);
    };

    if(!higgs_candidates.empty()) return *std::min_element(higgs_candidates.begin(), higgs_candidates.end(), Comparitor);

    return boost::optional<size_t>();
}

bool SignalObjectSelector::PassLeptonVetoSelection(const ntuple::Event& event)
{
    for(unsigned n = 0; n < event.other_lepton_p4.size(); ++n){
        if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::e) {
            const DiscriminatorIdResults eleId_iso(event.other_lepton_eleId_iso.at(n));
            const DiscriminatorIdResults eleId_noIso(event.other_lepton_eleId_noIso.at(n));
            const float iso = event.other_lepton_iso.at(n);
            if(eleId_iso.Passed(DiscriminatorWP::Medium)
                    || (eleId_noIso.Passed(DiscriminatorWP::Medium)
                        && iso < cuts::H_tautau_Run2::electronVeto::pfRelIso04))
                return false;
        }
        if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::mu) {
            DiscriminatorIdResults muonId(event.other_lepton_muonId.at(n));
            if((muonId.Passed(DiscriminatorWP::Medium) || muonId.Passed(DiscriminatorWP::Tight)) &&
                event.other_lepton_iso.at(n) < cuts::H_tautau_Run2::muonVeto::pfRelIso04) return false;
        }
    }
    return true;
}

bool SignalObjectSelector::PassMETfilters(const ntuple::Event& event, const Period period, bool is_Data)
{
    using Filter = ntuple::MetFilters::Filter;
    auto event_metFilters = ntuple::MetFilters(event.metFilters);
    if (!event_metFilters.Pass(Filter::PrimaryVertex)  || !event_metFilters.Pass(Filter::BeamHalo) ||
        !event_metFilters.Pass(Filter::HBHE_noise)  || !event_metFilters.Pass(Filter::HBHEiso_noise) ||
        !event_metFilters.Pass(Filter::ECAL_TP) || !event_metFilters.Pass(Filter::badMuon)) return false;
    if((period == Period::Run2017 || period == Period::Run2018) && !event_metFilters.Pass(Filter::ecalBadCalib))
        return false; //not in 2016
    if(is_Data && !event_metFilters.Pass(Filter::ee_badSC_noise)) return false;
    return true;
}

bool SignalObjectSelector::PassHTT_LeptonSelection(const LepCandidate& lepton, Channel channel, bool is_sync) const
{
    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_Run2::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_Run2::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_Run2::TauTau::tauID::pt}
      };

    auto e_id = GetTauVSeDiscriminator(channel);
    auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e)
        return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::H_tautau_Run2::MuTau::muonID::pt)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if(lepton->decayMode() == 5 || lepton->decayMode() == 6) return false;
    if(!is_sync && !lepton->Passed(e_id.first, e_id.second)) return false;
    if(!is_sync && !lepton->Passed(mu_id.first, mu_id.second)) return false;
    if(!lepton->Passed(GetTauVSjetDiscriminator().first, GetTauVSjetDiscriminator().second)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_LeptonSelection(const LepCandidate& lepton, Channel channel) const
{
    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_Run2::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_Run2::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_Run2::TauTau::tauID::pt}
      };

    auto e_id = GetTauVSeDiscriminator(channel);
    auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::H_tautau_Run2::MuTau::muonID::pt)) return false;
        if(!(lepton->iso() < cuts::H_tautau_Run2::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if(lepton->decayMode() == 5 || lepton->decayMode() == 6) return false;
    if(!lepton->Passed(e_id.first, e_id.second)) return false;
    if(!lepton->Passed(mu_id.first, mu_id.second)) return false;
    return true;
}

bool SignalObjectSelector::PassHH_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId,
                                                  bool is_sync) const
{
    static const std::map<Channel,double> pt_map_muon = {
        { Channel::MuTau, cuts::hh_bbtautau_Run2::MuTau::muonID::pt},
        { Channel::MuMu, cuts::hh_bbtautau_Run2::MuMu::muonID::pt},
    };

    static const std::map<Channel,double> pt_map_tau = {
        { Channel::ETau, cuts::hh_bbtautau_Run2::ETau::tauID::pt} ,
        { Channel::MuTau, cuts::hh_bbtautau_Run2::MuTau::tauID::pt},
        { Channel::TauTau, cuts::hh_bbtautau_Run2::TauTau::tauID::pt },
    };

    static const std::map<Channel,double> eta_map_tau = {
        { Channel::ETau, cuts::hh_bbtautau_Run2::ETau::tauID::eta },
        { Channel::MuTau, cuts::hh_bbtautau_Run2::MuTau::tauID::eta },
        { Channel::TauTau, cuts::hh_bbtautau_Run2::TauTau::tauID::eta },
    };

    if(lepton->leg_type() == LegType::e) {
        return lepton->passEleIsoId(DiscriminatorWP::Tight);
    }
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > pt_map_muon.at(channel))) return false;
        if(!lepton->passMuonId(DiscriminatorWP::Tight)) return false;
        if(legId == 1 && !(lepton->iso() < cuts::hh_bbtautau_Run2::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw exception("Leg Type Default Selection not supported");

    const auto e_id = GetTauVSeDiscriminator(channel);
    const auto mu_id = GetTauVSmuDiscriminator(channel);

    if(!(lepton.GetMomentum().pt() > pt_map_tau.at(channel))) return false;
    if(!(std::abs(lepton.GetMomentum().eta()) < eta_map_tau.at(channel))) return false;
    if(mode == SignalMode::HH && !lepton->PassedNewDecayMode()) return false;
    if(mode == SignalMode::HH_legacy && !lepton->PassedOldDecayMode()) return false;
    if((mode == SignalMode::HH && (lepton->decayMode() == 5 || lepton->decayMode() == 6))) return false;
    if(!(std::abs(lepton->dz()) < cuts::hh_bbtautau_Run2::TauTau::tauID::dz)) return false;
    if(std::abs(lepton->charge()) != 1) return false;
    if(!lepton->Passed(e_id.first, e_id.second)) return false;
    if(!lepton->Passed(mu_id.first, mu_id.second)) return false;
    const DiscriminatorWP first_leg_id = is_sync ? DiscriminatorWP::VVVLoose : GetTauVSjetDiscriminator().second;
    if(legId == 1 && !lepton->Passed(GetTauVSjetDiscriminator().first, first_leg_id)) return false;
    if(legId == 2 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                     GetTauVSjetSidebandWPRange().first)) return false;

    return true;
}

SignalObjectSelector::SelectedSignalJets::SelectedSignalJets() :
    bjet_pair(ntuple::LegPair::Undefined), vbf_pair(ntuple::LegPair::Undefined), n_bjets(0)
{
}

bool SignalObjectSelector::SelectedSignalJets::HasBjetPair() const { return bjet_pair.IsDefined(); }
bool SignalObjectSelector::SelectedSignalJets::HasVBFPair() const { return vbf_pair.IsDefined(); }
bool SignalObjectSelector::SelectedSignalJets::isSelectedBjet(size_t n) const { return bjet_pair.Contains(n); }
bool SignalObjectSelector::SelectedSignalJets::isSelectedVBFjet(size_t n) const { return vbf_pair.Contains(n); }

SignalObjectSelector::JetInfoCollection SignalObjectSelector::CreateJetInfos(const EventCandidate& event_candidate,
                                                                             const BTagger& bTagger,
                                                                             bool apply_jet_up_id)
{
    return CreateJetInfos(event_candidate, bTagger, apply_jet_up_id, boost::optional<size_t>(), SelectedSignalJets());
}

SignalObjectSelector::JetInfoCollection SignalObjectSelector::CreateJetInfos(
        const EventCandidate& event_candidate, const BTagger& bTagger, bool apply_jet_up_id,
        const boost::optional<size_t>& selected_htt_index, const SelectedSignalJets& selected_signal_jets)
{
    static constexpr double DeltaR_Lep_Jet = cuts::hh_bbtautau_Run2::DeltaR_Lep_Jet;

    JetInfoCollection jet_info_vector;
    const auto& event = event_candidate.GetEvent();
    std::vector<LorentzVectorM> lep_p4;
    if(selected_htt_index.is_initialized()) {
        if(*selected_htt_index >= event.first_daughter_indexes.size())
            throw exception("selected_htt_index is out of range to index first_daughter_indexes branch.");
        const size_t first_leg_id = event.first_daughter_indexes.at(*selected_htt_index);
        if(first_leg_id >= event_candidate.GetLeptons().size())
            throw exception("first_leg_id is out of range.");
        const auto& first_leg = event_candidate.GetLeptons().at(first_leg_id).GetMomentum();
        lep_p4.emplace_back(first_leg);
        if(*selected_htt_index >= event.second_daughter_indexes.size())
            throw exception("selected_htt_index is out of range to index second_daughter_indexes branch.");
        const size_t second_leg_id = event.second_daughter_indexes.at(*selected_htt_index);
        if(second_leg_id >= event_candidate.GetLeptons().size())
            throw exception("second_leg_id is out of range.");
        const auto& second_leg = event_candidate.GetLeptons().at(second_leg_id).GetMomentum();
        lep_p4.emplace_back(second_leg);
    }

    for(size_t n = 0; n < event_candidate.GetJets().size(); ++n) {
        const auto& jet_p4 = event_candidate.GetJets().at(n).GetMomentum();
        bool pass_dR_leptons = true;
        for(size_t lep_index = 0; lep_index < lep_p4.size() && pass_dR_leptons; ++lep_index)
            pass_dR_leptons = ROOT::Math::VectorUtil::DeltaR(lep_p4.at(lep_index), jet_p4) > DeltaR_Lep_Jet;
        if(!pass_dR_leptons) continue;
        if(selected_signal_jets.isSelectedBjet(n)) continue;
        if(selected_signal_jets.isSelectedVBFjet(n)) continue;
        DiscriminatorIdResults jet_pu_id(event_candidate.GetJets().at(n)->GetPuId());
        if(!PassEcalNoiceVeto(jet_p4, event_candidate.GetPeriod(), jet_pu_id)) continue;
        if(apply_jet_up_id && jet_p4.pt() < 50 && !jet_pu_id.Passed(DiscriminatorWP::Loose)) continue;
        const double tag = bTagger.BTag(*event_candidate.GetJets().at(n), false);
        jet_info_vector.emplace_back(jet_p4, n, tag);
    }
    return jet_info_vector;
}

SignalObjectSelector::SelectedSignalJets SignalObjectSelector::SelectSignalJets(
        const EventCandidate& event_candidate, size_t selected_htt_index, const BTagger& bTagger,
        DiscriminatorWP btag_wp)
{
    static constexpr bool apply_jet_up_id = true;
    const auto bjet_pt_cut = Cut1D_Bound::L(bTagger.PtCut());
    const auto bjet_eta_cut = Cut1D_Bound::AbsU(bTagger.EtaCut());

    const BTagger ptTagger(event_candidate.GetPeriod(), BTaggerKind::Pt);

    SelectedSignalJets selected_signal_jets;

    auto _CreateJetInfos = [&](const BTagger& tagger) {
        return CreateJetInfos(event_candidate, tagger, apply_jet_up_id, selected_htt_index, selected_signal_jets);
    };

    auto jet_info_vector = _CreateJetInfos(bTagger);
    auto bjets_ordered = jet_ordering::OrderJets(jet_info_vector, true, bjet_pt_cut, bjet_eta_cut);
    selected_signal_jets.n_bjets = bjets_ordered.size();
    if(bjets_ordered.size() >= 1){
        selected_signal_jets.bjet_pair.first = bjets_ordered.at(0).index;
    }

    if(bjets_ordered.size() >= 2){
        if(bTagger.GetTagger() == BTaggerKind::HHbtag
                || bTagger.Pass(*event_candidate.GetJets().at(bjets_ordered.at(1).index), btag_wp)) {
            selected_signal_jets.bjet_pair.second = bjets_ordered.at(1).index;
        }
    }
    auto jet_info_vector_vbf = _CreateJetInfos(ptTagger);
    auto vbf_jets_ordered = jet_ordering::OrderJets(jet_info_vector_vbf, true,
                                                    Cut1D_Bound::L(cuts::hh_bbtautau_Run2::jetID::vbf_pt),
                                                    Cut1D_Bound::AbsU(cuts::hh_bbtautau_Run2::jetID::vbf_eta));

    double max_mjj = -std::numeric_limits<double>::infinity();
    for(size_t n = 0; n < vbf_jets_ordered.size(); ++n) {
        const auto& jet_1 = vbf_jets_ordered.at(n);
        for(size_t h = n+1; h < vbf_jets_ordered.size(); ++h) {
            const auto& jet_2 = vbf_jets_ordered.at(h);
            const auto jet_12 = jet_1.p4 + jet_2.p4;
            if(jet_12.M() > max_mjj){
                max_mjj = jet_12.M();
                selected_signal_jets.vbf_pair = std::make_pair(vbf_jets_ordered.at(n).index,
                                                               vbf_jets_ordered.at(h).index);
            }
        }
    }

    if(selected_signal_jets.HasBjetPair()) return selected_signal_jets;

    auto jet_info_vector_new = _CreateJetInfos(bTagger);
    auto new_bjets_ordered = jet_ordering::OrderJets(jet_info_vector_new, true, bjet_pt_cut, bjet_eta_cut);
    if(new_bjets_ordered.size() >= 1)
        selected_signal_jets.bjet_pair.second = new_bjets_ordered.at(0).index;
    else{
        selected_signal_jets.vbf_pair = ntuple::LegPair::Undefined;
        if (bjets_ordered.size() >= 2)
            selected_signal_jets.bjet_pair.second = bjets_ordered.at(1).index;
    }

    return selected_signal_jets;
}

const FatJetCandidate* SignalObjectSelector::SelectFatJet(const EventCandidate& event_candidate,
                                                          const SelectedSignalJets& selected_signal_jets)
{
    using namespace cuts::hh_bbtautau_Run2;

    if(!selected_signal_jets.HasBjetPair()) return nullptr;
    for(const auto& fatJet : event_candidate.GetFatJets()) {
        if(fatJet->m(ntuple::TupleFatJet::MassType::SoftDrop) < fatJetID::mass) continue;
        if(fatJet->subJets().size() < 2) continue;
        auto subJets = fatJet->subJets();
        std::sort(subJets.begin(), subJets.end(), [](const ntuple::TupleSubJet& j1, const ntuple::TupleSubJet& j2) {
            return j1.p4().Pt() > j2.p4().Pt();
        });
        std::vector<double> deltaR;
        for(size_t n = 0; n < 2; ++n) {
            for(size_t k = 0; k < 2; ++k) {
                const auto& b_cand = event_candidate.GetJets().at(selected_signal_jets.bjet_pair.Get(k+1));
                const auto dR = ROOT::Math::VectorUtil::DeltaR(subJets.at(n).p4(), b_cand.GetMomentum());
                deltaR.push_back(dR);
            }
        }
        if((deltaR.at(0) < fatJetID::deltaR_subjet && deltaR.at(3) < fatJetID::deltaR_subjet)
                || (deltaR.at(1) < fatJetID::deltaR_subjet && deltaR.at(2) < fatJetID::deltaR_subjet))
            return &fatJet;
    }
    return nullptr;

}

bool SignalObjectSelector::PassEcalNoiceVetoImpl(const LorentzVector& jet_p4, Period period,
                                                 DiscriminatorIdResults jet_pu_id)
{
    if(period != Period::Run2017)
        return true;

    const double abs_eta = std::abs(jet_p4.eta());
    return !(jet_p4.pt() < cuts::hh_bbtautau_Run2::jetID::max_pt_veto
            && abs_eta > cuts::hh_bbtautau_Run2::jetID::eta_low_veto
            && abs_eta < cuts::hh_bbtautau_Run2::jetID::eta_high_veto
            && !jet_pu_id.Passed(DiscriminatorWP::Loose));
}

} // namespace analysis
