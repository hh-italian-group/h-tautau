/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"

namespace analysis {

SignalObjectSelector::SignalObjectSelector(SignalMode _mode) : mode(_mode)
{
    if(mode == SignalMode::HTT || mode == SignalMode::HTT_sync){
	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
    else if(mode == SignalMode::TauPOG_default){
    	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
    else if(mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full){
	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byDeepTau2017v2p1VSjet;
    }
    else if(mode == SignalMode::HH_legacy){
    	DR2_leptons = std::pow(cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
    else if(mode == SignalMode::HH){
    	 DR2_leptons = std::pow(cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects, 2);
	     discriminator = TauIdDiscriminator::byDeepTau2017v2p1VSjet;
    }

    else if(mode == SignalMode::Skimmer || mode == SignalMode::TauPOG_Skimmer){
        double DR = std::min(cuts::H_tautau_2016::DeltaR_betweenSignalObjects,cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects);
        DR2_leptons = std::pow(DR, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
    else
        throw analysis::exception("Signal Mode for SignalObjectSelector constructor not supported");
}

bool SignalObjectSelector::PassLeptonSelection(const LepCandidate& lepton, Channel channel, const size_t legId) const
{
    if(mode == SignalMode::HTT || mode == SignalMode::HTT_sync)
        return PassHTT_LeptonSelection(lepton,channel,mode == SignalMode::HTT_sync);
    if(mode == SignalMode::TauPOG_default || mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full)
	return PassTauPOG_LeptonSelection(lepton,channel);
    if(mode == SignalMode::HH_legacy)
        return PassHH_legacy_LeptonSelection(lepton,channel,legId);
    if(mode == SignalMode::HH)
        return PassHH_LeptonSelection(lepton,channel,legId);
    if(mode == SignalMode::Skimmer)
        return PassSkimmer_LeptonSelection(lepton);
    if(mode == SignalMode::TauPOG_Skimmer)
        return PassTauPOG_Skimmer_LeptonSelection(lepton);
    throw analysis::exception("Signal Mode for SignalObjectSelector class not supported");
}

boost::optional<size_t> SignalObjectSelector::GetHiggsCandidateIndex(EventCandidate& event_candidate) const
{
    const ntuple::Event& event = event_candidate.GetEvent();
    // std::vector<ntuple::TupleLepton> lepton_candidates;
    // for(size_t n = 0; n < event.lep_p4.size(); ++n)
    //     lepton_candidates.emplace_back(event, n);
    std::vector<size_t> higgs_candidates;
    for(size_t n = 0; n < event.first_daughter_indexes.size(); ++n){
        if(event.first_daughter_indexes.at(n) >=event_candidate.GetLeptons().size() ||
            event.second_daughter_indexes.at(n) >=event_candidate.GetLeptons().size())
            throw analysis::exception("Daughter indexes greater than lepton size.");
        auto& first_leg = event_candidate.GetLeptons().at(event.first_daughter_indexes.at(n));
        auto& second_leg = event_candidate.GetLeptons().at(event.second_daughter_indexes.at(n));
        if(ROOT::Math::VectorUtil::DeltaR2(first_leg.GetMomentum(), second_leg.GetMomentum()) <= DR2_leptons) continue;
        const Channel channel = static_cast<Channel>(event.channelId);
        if(!PassLeptonSelection(first_leg,channel,1)) continue;
        if(!PassLeptonSelection(second_leg,channel,2)) continue;
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
                const int iso_cmp = h1_leg->CompareIsolations(*h2_leg, discriminator);
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
            const Channel channel = static_cast<Channel>(event.channelId);
            if(channel == Channel::MuTau){ //this is a temporary fix because of an error in the production code for vetoMuons
                if(static_cast<LegType>(event.lep_type.at(0)) != LegType::mu)
                    throw analysis::exception("First leg type is not a muon in mutau channel.");
                if(ROOT::Math::VectorUtil::DeltaR2(event.lep_gen_p4.at(0), event.other_lepton_p4.at(n)) <= 0.1) continue;
            }
            analysis::DiscriminatorIdResults muonId(event.other_lepton_muonId.at(n));
            if(muonId.Passed(DiscriminatorWP::Medium)) return false;
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
    //againstElectron first, againstMuon second
    static const std::map<Channel,std::pair<DiscriminatorWP,DiscriminatorWP>> againstDiscriminators =
        { {Channel::ETau,{DiscriminatorWP::Tight,DiscriminatorWP::Loose}} ,
          {Channel::MuTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Tight}},
          {Channel::TauTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Loose}}
      };

    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_2016::TauTau::tauID::pt}
      };

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
    if(!is_sync && !lepton->Passed(TauIdDiscriminator::againstElectronMVA6,againstDiscriminators.at(channel).first)) return false;
    if(!is_sync && !lepton->Passed(TauIdDiscriminator::againstMuon3,againstDiscriminators.at(channel).second)) return false;
    if(!lepton->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_LeptonSelection(const LepCandidate& lepton, Channel channel) const
{
    //againstElectron first, againstMuon second
    static const std::map<Channel,std::pair<DiscriminatorWP,DiscriminatorWP>> againstDiscriminators =
        { {Channel::ETau,{DiscriminatorWP::Tight,DiscriminatorWP::Loose}} ,
          {Channel::MuTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Tight}},
          {Channel::TauTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Loose}}
      };

    // WP for deepTau vs E (first) and deepTau vs Mu (second)
    static const std::map<Channel,std::pair<DiscriminatorWP,DiscriminatorWP>> deepTauDiscriminators =
        { {Channel::ETau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VVVLoose}} ,
          {Channel::MuTau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VVVLoose}},
          {Channel::TauTau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VVVLoose}}
      };

    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_2016::TauTau::tauID::pt}
      };

    if(lepton->leg_type() == LegType::e) return true;
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::H_tautau_2017::MuTau::muonID::pt)) return false; //to be back
        //if(!(lepton.iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if((mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full) && (lepton->decayMode() == 5 || lepton->decayMode() == 6)) return false;
    if((mode == SignalMode::TauPOG_default) && !(lepton->PassedOldDecayMode())) return false;
    TauIdDiscriminator eleDiscriminator = mode == SignalMode::TauPOG_deepTauVsJet_full ? TauIdDiscriminator::byDeepTau2017v2p1VSe : TauIdDiscriminator::againstElectronMVA6;
    TauIdDiscriminator muonDiscriminator = mode == SignalMode::TauPOG_deepTauVsJet_full ? TauIdDiscriminator::byDeepTau2017v2p1VSmu : TauIdDiscriminator::againstMuon3;
    DiscriminatorWP eleWP = mode == SignalMode::TauPOG_deepTauVsJet_full ? deepTauDiscriminators.at(channel).first : againstDiscriminators.at(channel).first;
    DiscriminatorWP muonWP = mode == SignalMode::TauPOG_deepTauVsJet_full ? deepTauDiscriminators.at(channel).second : againstDiscriminators.at(channel).second;
    if(!lepton->Passed(eleDiscriminator,eleWP)) return false;
    if(!lepton->Passed(muonDiscriminator,muonWP)) return false;
    return true;
}

bool SignalObjectSelector::PassHH_legacy_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId) const
{
    //againstElectron first, againstMuon second
    static const std::map<Channel,std::pair<DiscriminatorWP,DiscriminatorWP>> againstDiscriminators =
        { {Channel::ETau,{DiscriminatorWP::Tight,DiscriminatorWP::Loose}} ,
          {Channel::MuTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Tight}},
          {Channel::TauTau,{DiscriminatorWP::VLoose,DiscriminatorWP::Loose}}
      };

    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::hh_bbtautau_2017::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::hh_bbtautau_2017::MuTau::tauID::pt},
          {Channel::TauTau, cuts::hh_bbtautau_2017::TauTau::tauID::pt}
      };

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
    if(!lepton->Passed(TauIdDiscriminator::againstElectronMVA6,againstDiscriminators.at(channel).first)) return false;
    if(!lepton->Passed(TauIdDiscriminator::againstMuon3,againstDiscriminators.at(channel).second)) return false;
    if(legId == 1 && !lepton->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::Medium)) return false;
    if(legId == 2 && !lepton->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
    return true;
}

//to be changed
bool SignalObjectSelector::PassHH_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId) const
{

    // WP for deepTau vs E (first) and deepTau vs Mu (second)
    static const std::map<Channel,std::pair<DiscriminatorWP,DiscriminatorWP>> deepTauDiscriminators =
        { {Channel::ETau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VLoose}} ,
          {Channel::MuTau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VLoose}},
          {Channel::TauTau,{DiscriminatorWP::VVVLoose,DiscriminatorWP::VLoose}}
      };

    static const std::map<Channel,double> pt_map =
        { {Channel::ETau, cuts::H_tautau_2016::ETau::tauID::pt} ,
          {Channel::MuTau, cuts::H_tautau_2016::MuTau::tauID::pt},
          {Channel::TauTau, cuts::H_tautau_2016::TauTau::tauID::pt}
      };

    if(lepton->leg_type() == LegType::e) {
        if(!lepton->passConversionVeto()) return false;
        if(!lepton->passEleIso(DiscriminatorWP::Medium)) return false;
        return true;
    }
    if(lepton->leg_type() == LegType::mu) {
        if(!(lepton.GetMomentum().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        if(legId == 1 && !(lepton->iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.GetMomentum().pt() > pt_map.at(channel))) return false;
    if((mode == SignalMode::HH) && !(lepton->PassedNewDecayMode())) return false;
    if((mode == SignalMode::HH && (lepton->decayMode() == 5 || lepton->decayMode() == 6))) return false;
    if(!lepton->Passed(TauIdDiscriminator::byDeepTau2017v2p1VSe,deepTauDiscriminators.at(channel).first)) return false;
    if(!lepton->Passed(TauIdDiscriminator::byDeepTau2017v2p1VSmu,deepTauDiscriminators.at(channel).second)) return false;
    if(legId == 1 && !lepton->Passed(TauIdDiscriminator::byDeepTau2017v2p1VSjet,DiscriminatorWP::Medium)) return false;
    if(legId == 2 && !lepton->Passed(TauIdDiscriminator::byDeepTau2017v2p1VSjet,DiscriminatorWP::VVVLoose)) return false;
    return true;
}

bool SignalObjectSelector::PassSkimmer_LeptonSelection(const LepCandidate& lepton) const
{
    if(lepton->leg_type() == LegType::e || lepton->leg_type() == LegType::mu) return true;
    if(!(lepton->leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton->PassedOldDecayMode())) return false;
    if(!lepton->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
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
                                                                  size_t selected_higgs_index)
{
    BTagger bTagger(period, jet_ordering);
    const double bjet_pt_cut = bTagger.PtCut();
    const double bjet_eta_cut = bTagger.EtaCut();

    SelectedSignalJets selected_signal_jets;
    const ntuple::Event& event = event_candidate.GetEvent();

    const auto CreateJetInfo = [&](bool useBTag) -> auto {
        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        for(size_t n = 0; n < event_candidate.GetJets().size(); ++n) {
            const size_t first_leg_id = event.first_daughter_indexes.at(selected_higgs_index);
            const auto& first_leg = event_candidate.GetLeptons().at(first_leg_id).GetMomentum();
            if(ROOT::Math::VectorUtil::DeltaR(first_leg, event_candidate.GetJets().at(n).GetMomentum()) <= cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
            const size_t second_leg_id = event.second_daughter_indexes.at(selected_higgs_index);
            const auto& second_leg = event_candidate.GetLeptons().at(second_leg_id).GetMomentum();
            if(ROOT::Math::VectorUtil::DeltaR(second_leg, event_candidate.GetJets().at(n).GetMomentum()) <= cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
            if(selected_signal_jets.isSelectedBjet(n)) continue;
            if(selected_signal_jets.isSelectedVBFjet(n)) continue;
            analysis::DiscriminatorIdResults jet_pu_id(event_candidate.GetJets().at(n)->GetPuId());
            if(!PassEcalNoiceVetoJets(event_candidate.GetJets().at(n).GetMomentum(), period, jet_pu_id)) continue;
            if(!jet_pu_id.Passed(analysis::DiscriminatorWP::Loose)) continue;
//            if(useBTag && (event.jets_pu_id.at(n) & (1 << 2)) == 0) continue;

            const double tag = useBTag ? bTagger.BTag(event,n) : event_candidate.GetJets().at(n).GetMomentum().Pt();
            jet_info_vector.emplace_back(event_candidate.GetJets().at(n).GetMomentum(),n,tag);
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
        selected_signal_jets.selectedVBFjetPair = ntuple::UndefinedLegPair();
        if (bjets_ordered.size() >= 2)
            selected_signal_jets.selectedBjetPair.second = bjets_ordered.at(1).index;
    }
    return selected_signal_jets;
}

} // namespace analysis
