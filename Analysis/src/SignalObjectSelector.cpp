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
	discriminator = TauIdDiscriminator::byDeepTau2017v1VSjet;
    }
    else if(mode == SignalMode::TauPOG_dpfTau){
    	DR2_leptons = std::pow(cuts::H_tautau_2016::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byDpfTau2016v0VSall;
    }
    else if(mode == SignalMode::HH){
    	DR2_leptons = std::pow(cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
        
    else if(mode == SignalMode::Skimmer || mode == SignalMode::TauPOG_Skimmer){
        double DR = std::min(cuts::H_tautau_2016::DeltaR_betweenSignalObjects,cuts::hh_bbtautau_2017::DeltaR_betweenSignalObjects);
        DR2_leptons = std::pow(DR, 2);
	discriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
    }
    else
        throw analysis::exception("Signal Mode for SignalObjectSelector constructor not supported");
}

bool SignalObjectSelector::PassLeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const
{
    if(mode == SignalMode::HTT || mode == SignalMode::HTT_sync)
        return PassHTT_LeptonSelection(lepton,channel,mode == SignalMode::HTT_sync);
    if(mode == SignalMode::TauPOG_default || mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full || 
	mode == SignalMode::TauPOG_dpfTau)
	return PassTauPOG_LeptonSelection(lepton,channel);
    if(mode == SignalMode::HH)
        return PassHH_LeptonSelection(lepton,channel);
    if(mode == SignalMode::Skimmer)
        return PassSkimmer_LeptonSelection(lepton);
    if(mode == SignalMode::TauPOG_Skimmer)
        return PassTauPOG_Skimmer_LeptonSelection(lepton);
    throw analysis::exception("Signal Mode for SignalObjectSelector class not supported");
}

boost::optional<size_t> SignalObjectSelector::GetHiggsCandidateIndex(const ntuple::Event& event) const
{
    std::vector<ntuple::TupleLepton> lepton_candidates;
    for(size_t n = 0; n < event.lep_p4.size(); ++n)
        lepton_candidates.emplace_back(event, n);

    std::vector<size_t> higgs_candidates;
    for(size_t n = 0; n < event.first_daughter_indexes.size(); ++n){
        const auto& first_leg = lepton_candidates.at(event.first_daughter_indexes.at(n));
        const auto& second_leg = lepton_candidates.at(event.second_daughter_indexes.at(n));
        if(ROOT::Math::VectorUtil::DeltaR2(first_leg.p4(), second_leg.p4()) <= DR2_leptons) continue;
        const Channel channel = static_cast<Channel>(event.channelId);
        if(!PassLeptonSelection(first_leg,channel)) continue;
        if(!PassLeptonSelection(second_leg,channel)) continue;
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
                const auto& h1_leg = lepton_candidates.at(h1_leg_id);
                const auto& h2_leg = lepton_candidates.at(h2_leg_id);
                const int iso_cmp = h1_leg.CompareIsolations(h2_leg, discriminator);
                if(iso_cmp != 0) return iso_cmp == 1;
                if(h1_leg.p4().pt() != h2_leg.p4().pt())
                    return h1_leg.p4().pt() > h2_leg.p4().pt();
            }
        }
        if(are_identical) return false;
        throw analysis::exception("not found a good criteria for best tau pair for %1%") % EventIdentifier(event);
    };

    if(!higgs_candidates.empty()) return *std::min_element(higgs_candidates.begin(), higgs_candidates.end(), Comparitor);
    return boost::optional<size_t>();
}

bool SignalObjectSelector::PassHTT_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel, bool is_sync) const
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

    if(lepton.leg_type() == LegType::e) return true;
    if(lepton.leg_type() == LegType::mu) {
        if(!(lepton.p4().pt() > cuts::H_tautau_2017::MuTau::muonID::pt)) return false; //to be back
	//if(!(lepton.p4().pt() > 30)) return false; //compatible with tau ID SF 
        //if(!is_sync && !(lepton.iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton.leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.p4().pt() > pt_map.at(channel))) return false;
    //if(lepton.decayMode() == 5 || lepton.decayMode() == 6 || lepton.decayMode() == 11) return false;
    if(!(lepton.PassedOldDecayMode())) return false;
    if(!is_sync && !lepton.Passed(TauIdDiscriminator::againstElectronMVA6,againstDiscriminators.at(channel).first)) return false;
    if(!is_sync && !lepton.Passed(TauIdDiscriminator::againstMuon3,againstDiscriminators.at(channel).second)) return false;
    if(!lepton.Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const
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

    if(lepton.leg_type() == LegType::e) return true;
    if(lepton.leg_type() == LegType::mu) {
        if(!(lepton.p4().pt() > cuts::H_tautau_2017::MuTau::muonID::pt)) return false; //to be back
        //if(!(lepton.iso() < cuts::H_tautau_2016::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton.leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.p4().pt() > pt_map.at(channel))) return false;
    if((mode == SignalMode::TauPOG_deepTauVsJet || mode == SignalMode::TauPOG_deepTauVsJet_full) && (lepton.decayMode() == 5 || lepton.decayMode() == 6)) return false;
    if((mode == SignalMode::TauPOG_default || mode == SignalMode::TauPOG_dpfTau) && !(lepton.PassedOldDecayMode())) return false;
    TauIdDiscriminator eleDiscriminator = mode == SignalMode::TauPOG_deepTauVsJet_full ? TauIdDiscriminator::byDeepTau2017v1VSe : 												 TauIdDiscriminator::againstElectronMVA6;
    TauIdDiscriminator muonDiscriminator = mode == SignalMode::TauPOG_deepTauVsJet_full ? TauIdDiscriminator::byDeepTau2017v1VSmu : 												  TauIdDiscriminator::againstMuon3;
    DiscriminatorWP eleWP = mode == SignalMode::TauPOG_deepTauVsJet_full ? deepTauDiscriminators.at(channel).first : 										   againstDiscriminators.at(channel).first;
    DiscriminatorWP muonWP = mode == SignalMode::TauPOG_deepTauVsJet_full ? deepTauDiscriminators.at(channel).second : 										    againstDiscriminators.at(channel).second;
    if(!lepton.Passed(eleDiscriminator,eleWP)) return false;
    if(!lepton.Passed(muonDiscriminator,muonWP)) return false;
    return true;
}

bool SignalObjectSelector::PassHH_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const
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

    if(lepton.leg_type() == LegType::e) return true;
    if(lepton.leg_type() == LegType::mu) {
        if(!(lepton.p4().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        if(!(lepton.iso() < cuts::hh_bbtautau_2017::MuTau::muonID::pfRelIso04)) return false;
        return true;
    }
    if(!(lepton.leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.p4().pt() > pt_map.at(channel))) return false;
    //if(lepton.decayMode() == 5 || lepton.decayMode() == 6 || lepton.decayMode() == 11) return false;
    if(!(lepton.PassedOldDecayMode())) return false;
    if(!lepton.Passed(TauIdDiscriminator::againstElectronMVA6,againstDiscriminators.at(channel).first)) return false;
    if(!lepton.Passed(TauIdDiscriminator::againstMuon3,againstDiscriminators.at(channel).second)) return false;
    if(!lepton.Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
    return true;
}

bool SignalObjectSelector::PassSkimmer_LeptonSelection(const ntuple::TupleLepton& lepton) const
{
    if(lepton.leg_type() == LegType::e || lepton.leg_type() == LegType::mu) return true;
    if(!(lepton.leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    if(!(lepton.PassedOldDecayMode())) return false;
    if(!lepton.Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,DiscriminatorWP::VVLoose)) return false;
    return true;
}

bool SignalObjectSelector::PassTauPOG_Skimmer_LeptonSelection(const ntuple::TupleLepton& lepton) const
{
    if(lepton.leg_type() == LegType::e) return true;
    if(lepton.leg_type() == LegType::mu) {
        if(!(lepton.p4().pt() > cuts::hh_bbtautau_2017::MuTau::muonID::pt)) return false;
        return true;
    }
    if(!(lepton.leg_type() == LegType::tau)) throw analysis::exception("Leg Type Default Selection not supported");
    return true;
}

} // namespace analysis
