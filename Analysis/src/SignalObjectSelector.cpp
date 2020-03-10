/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"

namespace analysis {

SignalObjectSelector::SignalObjectSelector(SignalMode _mode) : mode(_mode)
{
    if(mode == SignalMode::HTT || mode == SignalMode::HTT_sync)
	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
    else if(mode == SignalMode::TauPOG_default)
    	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
    else if(mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full)
	   DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
    else if(mode == SignalMode::HH_legacy)
    	DR2_leptons = std::pow(cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects, 2);
    else if(mode == SignalMode::HH)
    	 DR2_leptons = std::pow(cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects, 2);
    else if(mode == SignalMode::Skimmer || mode == SignalMode::TauPOG_Skimmer){
        double DR = std::min(cuts::H_tautau_2016::DeltaR_betweenSignalObjects,cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects);
        DR2_leptons = std::pow(DR, 2);
    }
    else
        throw analysis::exception("Signal Mode for SignalObjectSelector constructor not supported");
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSjetDiscriminator() const
{
    static const std::map<SignalMode, std::pair<TauIdDiscriminator, DiscriminatorWP>> tauVSjetDiscriminator = {
        { SignalMode::HTT, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::HTT_sync, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_default, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_deepTauVsJet, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_deepTauVsJet_full, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium } },
        { SignalMode::HH_legacy, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::HH, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium } },
        { SignalMode::Skimmer, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_Skimmer, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, DiscriminatorWP::Medium } }
    };

    if(!tauVSjetDiscriminator.count(mode))
        throw exception("Mode'%1%' isn't found at the tauVSjetDiscriminator map") % mode;
    return tauVSjetDiscriminator.at(mode);
}

std::pair<DiscriminatorWP, DiscriminatorWP> SignalObjectSelector::GetTauVSjetSidebandWPRange() const
{
    static const std::map<SignalMode, std::pair<DiscriminatorWP, DiscriminatorWP>> tauVSjetWps = {
        { SignalMode::HTT, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::HTT_sync, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_default, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_deepTauVsJet, { DiscriminatorWP::VVVLoose, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_deepTauVsJet_full, { DiscriminatorWP::VVVLoose, DiscriminatorWP::Medium } },
        { SignalMode::HH_legacy, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::HH, { DiscriminatorWP::VVVLoose, DiscriminatorWP::Medium } },
        { SignalMode::Skimmer, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } },
        { SignalMode::TauPOG_Skimmer, { DiscriminatorWP::VVLoose, DiscriminatorWP::Medium } }
    };

    if(!tauVSjetWps.count(mode))
        throw exception("Mode'%1%' isn't found at the tauVSjetWps map") % mode;
    return tauVSjetWps.at(mode);
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSeDiscriminator(analysis::Channel channel) const
{
    static const std::map<std::pair<Channel, SignalMode>, std::pair<TauIdDiscriminator, DiscriminatorWP>> tauVSeDiscriminators = {
        { {Channel::ETau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::ETau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HTT}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::ETau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VVVLoose} },
        { {Channel::TauTau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VVVLoose} },
        { {Channel::ETau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VVVLoose} },
        { {Channel::TauTau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VVVLoose} },
        { {Channel::ETau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} },
        { {Channel::TauTau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} },
        { {Channel::ETau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstElectronMVA6, DiscriminatorWP::VLoose} },
        { {Channel::ETau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VLoose} },
        { {Channel::TauTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVLoose} },
        { {Channel::MuMu, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::VVVLoose} }
        // { {Channel::ETau, SignalMode::Skimmer}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        // { {Channel::ETau, SignalMode::TauPOG_Skimmer}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
    };

    auto disc_wp = tauVSeDiscriminators.find({channel, mode});
    if(disc_wp == tauVSeDiscriminators.end())
        throw exception("Channel: %1% and mode: %2% not found in the tauVSeDiscriminators map") % channel % mode;
    return disc_wp->second;
}

std::pair<TauIdDiscriminator, DiscriminatorWP> SignalObjectSelector::GetTauVSmuDiscriminator(analysis::Channel channel) const
{
    static const std::map<std::pair<Channel, SignalMode>, std::pair<TauIdDiscriminator, DiscriminatorWP>> tauVSmuDiscriminators = {
        { {Channel::ETau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HTT_sync}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HTT}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::TauPOG_default}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::TauPOG_deepTauVsJet}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::TauPOG_deepTauVsJet_full}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::MuTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HH_legacy}, {TauIdDiscriminator::againstMuon3, DiscriminatorWP::Loose} },
        { {Channel::ETau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::MuTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::Tight} },
        { {Channel::TauTau, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::VLoose} },
        { {Channel::MuMu, SignalMode::HH}, {TauIdDiscriminator::byDeepTau2017v2p1VSmu, DiscriminatorWP::VLoose} }
        // { {Channel::ETau, SignalMode::Skimmer}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
        // { {Channel::ETau, SignalMode::TauPOG_Skimmer}, {TauIdDiscriminator::byDeepTau2017v2p1VSe, DiscriminatorWP::Tight} },
    };
    auto disc_wp = tauVSmuDiscriminators.find({channel, mode});
    if(disc_wp == tauVSmuDiscriminators.end())
        throw exception("Channel: %1% and mode: %2% not found in the tauVSmuDiscriminators map") % channel % mode;
    return disc_wp->second;
}

bool SignalObjectSelector::PassLeptonSelection(const LepCandidate& lepton, Channel channel, const size_t legId, bool is_sync) const
{
    if(mode == SignalMode::HTT || mode == SignalMode::HTT_sync)
        return PassHTT_LeptonSelection(lepton,channel,mode == SignalMode::HTT_sync);
    if(mode == SignalMode::TauPOG_default || mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full)
        return PassTauPOG_LeptonSelection(lepton,channel);
    if(mode == SignalMode::HH_legacy)
        return PassHH_legacy_LeptonSelection(lepton,channel,legId);
    if(mode == SignalMode::HH)
        return PassHH_LeptonSelection(lepton,channel,legId, is_sync);
    if(mode == SignalMode::Skimmer)
        return PassSkimmer_LeptonSelection(lepton);
    if(mode == SignalMode::TauPOG_Skimmer)
        return PassTauPOG_Skimmer_LeptonSelection(lepton);
    throw analysis::exception("Signal Mode for SignalObjectSelector class not supported");
}

boost::optional<size_t> SignalObjectSelector::GetHiggsCandidateIndex(EventCandidate& event_candidate, bool is_sync) const
{
    const ntuple::Event& event = event_candidate.GetEvent();
    std::vector<size_t> higgs_candidates;
    for(size_t n = 0; n < event.first_daughter_indexes.size(); ++n){
        if(event.first_daughter_indexes.at(n) >=event_candidate.GetLeptons().size() ||
            event.second_daughter_indexes.at(n) >=event_candidate.GetLeptons().size())
            throw analysis::exception("Daughter indexes greater than lepton size.");
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
        throw analysis::exception("not found a good criteria for best tau pair for %1%") % EventIdentifier(event);
    };

    if(!higgs_candidates.empty()) return *std::min_element(higgs_candidates.begin(), higgs_candidates.end(), Comparitor);

    return boost::optional<size_t>();
}

bool SignalObjectSelector::PassLeptonVetoSelection(const ntuple::Event& event) const
{
    for(unsigned n = 0; n < event.other_lepton_p4.size(); ++n){
        if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::e){
          analysis::DiscriminatorIdResults eleId_iso(event.other_lepton_eleId_iso.at(n));
          if(eleId_iso.Passed(DiscriminatorWP::Medium)) return false;
        }
        if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::mu){
            analysis::DiscriminatorIdResults muonId(event.other_lepton_muonId.at(n));
            if((muonId.Passed(DiscriminatorWP::Medium) || muonId.Passed(DiscriminatorWP::Tight)) &&
                event.other_lepton_iso.at(n) < 0.3) return false;
        }
    }
    return true;
}

bool SignalObjectSelector::PassMETfilters(const ntuple::Event& event, const analysis::Period period, bool is_Data) const
{
    using Filter = ntuple::MetFilters::Filter;
    auto event_metFilters = ntuple::MetFilters(event.metFilters);
    if (!event_metFilters.Pass(Filter::PrimaryVertex)  || !event_metFilters.Pass(Filter::BeamHalo) ||
        !event_metFilters.Pass(Filter::HBHE_noise)  || !event_metFilters.Pass(Filter::HBHEiso_noise) ||
        !event_metFilters.Pass(Filter::ECAL_TP) || !event_metFilters.Pass(Filter::badMuon)) return false;
    if((period == Period::Run2017 || period == Period::Run2018) && !event_metFilters.Pass(Filter::ecalBadCalib)) return false; //not in 2016
    if(is_Data && !event_metFilters.Pass(Filter::ee_badSC_noise)) return false;
    return true;
}

bool SignalObjectSelector::PassHTT_LeptonSelection(const LepCandidate& lepton, Channel channel, bool is_sync) const
{
    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_2016::TauTau::tauID::pt}
      };

    auto e_id = GetTauVSeDiscriminator(channel);
    auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::H_tautau_2017::MuTau::muonID::pt)) return false; //to be back
	//if(!(lepton.p4().pt() > 30)) return false; //compatible with tau ID SF
        //if(!is_sync && !(lepton.iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    //if(lepton.decayMode() == 5 || lepton.decayMode() == 6 || lepton.decayMode() == 11) return false;
    if(!(lepton->PassedOldDecayMode())) return false;
    if(!is_sync && !lepton->Passed(e_id.first, e_id.second)) return false;
    if(!is_sync && !lepton->Passed(mu_id.first, mu_id.second)) return false;
    if(!lepton->Passed(GetTauVSjetDiscriminator().first,GetTauVSjetDiscriminator().second)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_LeptonSelection(const LepCandidate& lepton, Channel channel) const
{
    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_2016::TauTau::tauID::pt}
      };

    auto e_id = GetTauVSeDiscriminator(channel);
    auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::H_tautau_2017::MuTau::muonID::pt)) return false; //to be back
        //if(!(lepton.iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if((mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full)
        && (lepton->decayMode() == 5 || lepton->decayMode() == 6)) return false;
    if((mode == SignalMode::TauPOG_default) && !(lepton->PassedOldDecayMode())) return false;
    if(!lepton->Passed(e_id.first,e_id.second)) return false;
    if(!lepton->Passed(mu_id.first,mu_id.second)) return false;
    return true;
}

bool SignalObjectSelector::PassHH_legacy_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId) const
{
    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::hh_bbtautau_2017::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::hh_bbtautau_2017::MuTau::tauID::pt},
          {Channel::TauTau, cuts::hh_bbtautau_2017::TauTau::tauID::pt}
      };

    auto e_id = GetTauVSeDiscriminator(channel);
    auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        if(legId == 1 && !(lepton->iso() < cuts::hh_bbtautau_2017::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    //if(lepton.decayMode() == 5 || lepton.decayMode() == 6 || lepton.decayMode() == 11) return false;
    if(!(lepton->PassedOldDecayMode())) return false;
    if(!lepton->Passed(e_id.first, e_id.second)) return false;
    if(!lepton->Passed(mu_id.first, mu_id.second)) return false;
    if(legId == 1 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                     GetTauVSjetDiscriminator().second)) return false;
    if(legId == 2 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                     GetTauVSjetSidebandWPRange().first)) return false;
    return true;
}

//to be changed
bool SignalObjectSelector::PassHH_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId, bool is_sync) const
{
    static const std::map<Channel,double> pt_map =
        { { Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          { Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          { Channel::TauTau, cuts::hh_bbtautau_2017::TauTau::tauID::pt },
      };

      static const std::map<Channel,double> eta_map =
          { { Channel::ETau, cuts::H_tautau_2016::ETau::tauID::eta },
            { Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::eta },
            { Channel::TauTau, cuts::hh_bbtautau_2017::TauTau::tauID::eta_sel },
        };

     auto e_id = GetTauVSeDiscriminator(channel);
     auto mu_id = GetTauVSmuDiscriminator(channel);

    if(lepton->leg_type() == LegType::e) {
        if(!lepton->passEleIso(DiscriminatorWP::Tight)) return false;
        return true;
    }
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        if(legId == 1 && !(lepton->iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        if(legId == 1 && !lepton->passMuonId(DiscriminatorWP::Tight)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if(!(std::abs(lepton.GetMomentum().eta()) < eta_map.at(channel))) return false;
    if((mode == SignalMode::HH) && !(lepton->PassedNewDecayMode())) return false;
    if((mode == SignalMode::HH && (lepton->decayMode() == 5 || lepton->decayMode() == 6))) return false;
    if(!lepton->Passed(e_id.first, e_id.second)) return false;
    if(!lepton->Passed(mu_id.first, mu_id.second)) return false;
    if(is_sync && legId == 1 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                                DiscriminatorWP::VVVLoose)) return false;
    if(!is_sync && legId == 1 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                                 GetTauVSjetDiscriminator().second)) return false;
    if(legId == 2 && !lepton->Passed(GetTauVSjetDiscriminator().first,
                                     GetTauVSjetSidebandWPRange().first)) return false;

    return true;
}

bool SignalObjectSelector::PassSkimmer_LeptonSelection(const LepCandidate& lepton) const
{
    if(lepton->leg_type() == LegType::e || lepton->leg_type() == LegType::mu) return true;
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton->PassedOldDecayMode())) return false;
    if(!lepton->Passed(GetTauVSjetDiscriminator().first,
                       GetTauVSjetSidebandWPRange().first)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_Skimmer_LeptonSelection(const LepCandidate& lepton) const
{
    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    return true;
}

SignalObjectSelector::SelectedSignalJets::SelectedSignalJets() : selectedBjetPair(ntuple::UndefinedLegPair()),
    selectedVBFjetPair(ntuple::UndefinedLegPair()), n_bjets(0) { }

bool SignalObjectSelector::SelectedSignalJets::HasBjetPair(size_t njets) const
{
    return selectedBjetPair.first < njets && selectedBjetPair.second < njets;
}

bool SignalObjectSelector::SelectedSignalJets::HasVBFPair(size_t njets) const
{
    return selectedVBFjetPair.first < njets && selectedVBFjetPair.second < njets;
}

bool SignalObjectSelector::SelectedSignalJets::isSelectedBjet(size_t n) const
{
    return selectedBjetPair.first == n || selectedBjetPair.second == n;
}

bool SignalObjectSelector::SelectedSignalJets::isSelectedVBFjet(size_t n) const
{
    return selectedVBFjetPair.first == n || selectedVBFjetPair.second == n;
}

SignalObjectSelector::SelectedSignalJets SignalObjectSelector::SelectSignalJets(EventCandidate& event_candidate,
                                                                  const analysis::Period& period,
                                                                  analysis::JetOrdering jet_ordering,
                                                                  size_t selected_higgs_index,
                                                                  UncertaintySource uncertainty_source,
                                                                  UncertaintyScale scale)
{
    BTagger bTagger(period, jet_ordering);
    bool base_ordering = jet_ordering != JetOrdering::HHJetTag;
    const double bjet_pt_cut = bTagger.PtCut();
    const double bjet_eta_cut = bTagger.EtaCut();

    SelectedSignalJets selected_signal_jets;
    const ntuple::Event& event = event_candidate.GetEvent();
    const auto CreateJetInfo = [&](bool useBTag) -> auto {
        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        for(size_t n = 0; n < event_candidate.GetJets().size(); ++n) {
            if(selected_higgs_index >= event.first_daughter_indexes.size())
                throw exception("selected_higgs_index is out of range to index first_daughter_indexes branch.");
            const size_t first_leg_id = event.first_daughter_indexes.at(selected_higgs_index);
            if(first_leg_id >= event_candidate.GetLeptons().size())
                throw exception("first_leg_id is out of range.");
            const auto& first_leg = event_candidate.GetLeptons().at(first_leg_id).GetMomentum();
            if(ROOT::Math::VectorUtil::DeltaR(first_leg, event_candidate.GetJets().at(n).GetMomentum()) <=
                    cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
            if(selected_higgs_index >= event.second_daughter_indexes.size())
                throw exception("selected_higgs_index is out of range to index second_daughter_indexes branch.");
            const size_t second_leg_id = event.second_daughter_indexes.at(selected_higgs_index);
            if(second_leg_id >= event_candidate.GetLeptons().size())
                throw exception("second_leg_id is out of range.");
            const auto& second_leg = event_candidate.GetLeptons().at(second_leg_id).GetMomentum();
            if(ROOT::Math::VectorUtil::DeltaR(second_leg, event_candidate.GetJets().at(n).GetMomentum()) <=
                    cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
            if(selected_signal_jets.isSelectedBjet(n)) continue;
            if(selected_signal_jets.isSelectedVBFjet(n)) continue;
            analysis::DiscriminatorIdResults jet_pu_id(event_candidate.GetJets().at(n)->GetPuId());
            if(!PassEcalNoiceVetoJets(event_candidate.GetJets().at(n).GetMomentum(), period, jet_pu_id)) continue;
            if(event_candidate.GetJets().at(n).GetMomentum().pt() < 50 && !jet_pu_id.Passed(analysis::DiscriminatorWP::Loose)) continue;
            const double tag = useBTag ? bTagger.BTag(event, n, uncertainty_source, scale, base_ordering) : event_candidate.GetJets().at(n).GetMomentum().Pt();
            jet_info_vector.emplace_back(event_candidate.GetJets().at(n).GetMomentum(),n,tag);
        }
        return jet_info_vector;
    };
    auto jet_info_vector = CreateJetInfo(true);
    auto bjets_ordered = jet_ordering::OrderJets(jet_info_vector, true, bjet_pt_cut, bjet_eta_cut);
    selected_signal_jets.n_bjets = bjets_ordered.size();
    if(bjets_ordered.size() >= 1){
        selected_signal_jets.selectedBjetPair.first = bjets_ordered.at(0).index;
    }

    if(bjets_ordered.size() >= 2){
        if(bTagger.Pass(event, bjets_ordered.at(1).index, uncertainty_source, scale)
                || jet_ordering == analysis::JetOrdering::HHJetTag) {
            selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
        }
    }
    auto jet_info_vector_vbf = CreateJetInfo(false);
    auto vbf_jets_ordered = jet_ordering::OrderJets(jet_info_vector_vbf, true,
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
    auto new_bjets_ordered = jet_ordering::OrderJets(jet_info_vector_new, true, bjet_pt_cut, bjet_eta_cut);
    if(new_bjets_ordered.size() >= 1)
        selected_signal_jets.selectedBjetPair.second = new_bjets_ordered.at(0).index;
    else{
        selected_signal_jets.selectedVBFjetPair = ntuple::UndefinedLegPair();
        if (bjets_ordered.size() >= 2)
            selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
    }

    return selected_signal_jets;
}

} // namespace analysis
