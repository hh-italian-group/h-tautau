/* Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/McCorrections/include/LeptonWeights.h"

#include <boost/algorithm/string/case_conv.hpp>
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {
namespace mc_corrections {
namespace detail {

LeptonScaleFactors::LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInputSingle,
                                       const std::string& triggerInputCross) :
        idIso(new SF()), triggerSingle(new SF())
{
    idIso->init_ScaleFactor(idIsoInput);
    triggerSingle->init_ScaleFactor(triggerInputSingle);
    if(!triggerInputCross.empty()) {
        triggerCross = std::make_shared<SF>();
        triggerCross->init_ScaleFactor(triggerInputCross);
    }
}

double LeptonScaleFactors::GetIdIsoSF(const LorentzVector& p4, UncertaintyScale unc_scale) const
{
    return idIso->get_ScaleFactor(p4.pt(), p4.eta())
           + static_cast<int>(unc_scale) * idIso->get_ScaleFactorError(p4.pt(), p4.eta());
}

double LeptonScaleFactors::GetTriggerEff(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale) const
{
    return GetEff(p4, isData, unc_scale, *triggerSingle);
}

double LeptonScaleFactors::GetTriggerEffCross(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale) const
{
    if(!triggerCross)
        throw exception("LeptonScaleFactors: cross-trigger is not set.");
    return GetEff(p4, isData, unc_scale, *triggerCross);
}

bool LeptonScaleFactors::HasCrossTriggers() const { return triggerCross.get() != nullptr; }

double LeptonScaleFactors::GetEff(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale,
                                  SF& eff_provider)
{
    double eff, err=0;
    if(isData) {
        eff = eff_provider.get_EfficiencyData(p4.pt(), p4.eta());
        err = eff_provider.get_EfficiencyDataError(p4.pt(), p4.eta());
    } else {
        eff = eff_provider.get_EfficiencyMC(p4.pt(), p4.eta());
        err = eff_provider.get_EfficiencyMCError(p4.pt(), p4.eta());
    }
    return eff + static_cast<int>(unc_scale) * err;
}

ElectronScaleFactorPOG::ElectronScaleFactorPOG(const std::string& idInput, const std::string& isoInput) :
    id_hist(LoadWeight(idInput,"EGamma_SF2D")), iso_hist(LoadWeight(isoInput,"EGamma_SF2D"))
{
}

double ElectronScaleFactorPOG::GetTriggerSF() const { return 0.991; }

ElectronScaleFactorPOG::HistPtr ElectronScaleFactorPOG::LoadWeight(const std::string& weight_file_name,
                                                                   const std::string& hist_name)
{
    auto file = root_ext::OpenRootFile(weight_file_name);
    return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
}

MuonScaleFactorPOG::MuonScaleFactorPOG(const std::string& idInput, const std::string& isoInput,
                                       const std::string& triggerInput) :
    id_hist(LoadWeight(idInput,"NUM_MediumID_DEN_genTracks_pt_abseta")),
    iso_hist(LoadWeight(isoInput,"NUM_TightRelIso_DEN_MediumID_pt_abseta")),
    trigger_hist(LoadWeight( triggerInput,"IsoMu27_PtEtaBins/pt_abseta_ratio"))
 {
 }

MuonScaleFactorPOG::HistPtr MuonScaleFactorPOG::LoadWeight(const std::string& weight_file_name,
                                                           const std::string& hist_name)
{
    auto file = root_ext::OpenRootFile(weight_file_name);
    return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
}

} // namespace detail

LeptonWeights::LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_SingletriggerInput,
                             const std::string& electron_CrossTriggerInput, const std::string& muon_idIsoInput,
                             const std::string& muon_SingletriggerInput, const std::string& muon_CrossTriggerInput,
                             const std::string& _tauTriggerInput, const std::string& _tauVBFTriggerInput,
                             Period _period, bool _is_dm_binned) :
    electronSF(electron_idIsoInput, electron_SingletriggerInput, electron_CrossTriggerInput),
    muonSF(muon_idIsoInput, muon_SingletriggerInput, muon_CrossTriggerInput), tauTriggerInput(_tauTriggerInput),
    tauVBFTriggerInput(_tauVBFTriggerInput), period(_period), is_dm_binned(_is_dm_binned)
{
}

TauIDSFTool& LeptonWeights::GetTauIdProvider(TauIdDiscriminator discr, DiscriminatorWP wp)
{
    static const std::map<Period ,std::string>  period_label = {
        { Period::Run2016, "2016Legacy" }, { Period::Run2017, "2017ReReco" }, { Period::Run2018, "2018ReReco" }
    };
     auto& provider = tau_sf_providers[discr][wp];
     if(!provider) {
         static const std::string s = "by";
         const auto i = ToString(discr).find(s);
         const std::string discr_name = ToString(discr).erase(i, s.length());
         if(!period_label.count(period))
             throw exception("LeptonWeights::GetTauIdProvider: Period %1% is not supported.") % period;
         provider = std::make_shared<TauIDSFTool>(period_label.at(period), discr_name, ToString(wp), is_dm_binned,
                                                  false);
     }
     return *provider;
}

const tau_trigger::SFProvider& LeptonWeights::GetTauTriggerSFProvider(Channel channel, DiscriminatorWP wp, bool is_vbf)
{
        static const std::map<Channel, std::string> channel_names = {
            { Channel::ETau, "etau" }, { Channel::MuTau, "mutau" }, { Channel::TauTau, "ditau" }
        };
        auto& provider = tau_trigger_sf_providers[std::make_pair(channel,is_vbf)][wp];
        if(!provider) {
            if(!channel_names.count(channel))
                throw exception("LeptonWeights::GetTauTriggerSFProvider: channel %1% is not supported.") % channel;
            if(!is_vbf)
                provider = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, channel_names.at(channel),
                                                                     ToString(wp));
            else
                provider = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "ditauvbf", ToString(wp));
        }
        return *provider;
}

double LeptonWeights::GetLegIdIsoWeight(LepCandidate leg, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                                        DiscriminatorWP VSjet_wp,
                                        UncertaintySource unc_source, UncertaintyScale unc_scale)
{
    const UncertaintyScale current_scale = ApplyIdUncertaintyScale(leg->decayMode(), leg.GetMomentum().pt(),
                                                                   leg.GetMomentum().eta(),leg->gen_match(),
                                                                   leg->leg_type(), unc_source)
                                         ? unc_scale : UncertaintyScale::Central;

    if(leg->leg_type() == LegType::e) {
        return electronSF.GetIdIsoSF(leg.GetMomentum(), current_scale);
    } else if(leg->leg_type() == LegType::mu) {
        return muonSF.GetIdIsoSF(leg.GetMomentum(), current_scale);
    } else if(leg->leg_type() == LegType::tau) {
        double tau_id_weight_vs_mu = 1;
        double tau_id_weight_vs_ele = 1;
        double tau_id_weight = 1;

        const std::string scale = current_scale == UncertaintyScale::Central ? "" : ToString(current_scale);
        const int gen_match = static_cast<int>(leg->gen_match());

        if(leg->gen_match() == GenLeptonMatch::Electron || leg->gen_match() == GenLeptonMatch::TauElectron) {
            const auto& tauIdWeightVsEle = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSe, VSe_wp);
            tau_id_weight_vs_ele = tauIdWeightVsEle.getSFvsEta(std::abs(leg.GetMomentum().eta()), gen_match, scale);
        } else if(leg->gen_match() == GenLeptonMatch::Muon || leg->gen_match() == GenLeptonMatch::TauMuon) {
            const auto& tauIdWeightVsMu = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSmu, VSmu_wp);
            tau_id_weight_vs_mu = tauIdWeightVsMu.getSFvsEta(std::abs(leg.GetMomentum().eta()), gen_match, scale);
        } else if(leg->gen_match() == GenLeptonMatch::Tau) {
            auto& tauIdWeight = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSjet, VSjet_wp);

            const auto getTauSF = [&](const std::string& scale_str) {
                return is_dm_binned ? tauIdWeight.getSFvsDM(leg.GetMomentum().pt(), leg->decayMode(),
                                                            gen_match, scale_str)
                                    : tauIdWeight.getSFvsPT(leg.GetMomentum().pt(), gen_match, scale_str);
            };

            tau_id_weight = getTauSF(scale);
            if(current_scale != UncertaintyScale::Central
                    && (VSe_wp < DiscriminatorWP::VLoose || VSmu_wp < DiscriminatorWP::Tight)) {
                const double tau_id_weight_central = getTauSF("");
                const double err = std::abs(tau_id_weight_central - tau_id_weight);
                const double rel_err = leg.GetMomentum().pt() < 100 ? 0.03 : 0.15;
                tau_id_weight = tau_id_weight_central + static_cast<int>(current_scale)
                                * (err + rel_err * tau_id_weight_central);
            }
        }
        return tau_id_weight * tau_id_weight_vs_mu * tau_id_weight_vs_ele;
    } else {
        throw exception("LeptonWeights::GetLegIdIsoWeight: Leg type not found.");
    }
}

double LeptonWeights::GetIdIsoWeight(EventInfo& eventInfo, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                                     DiscriminatorWP VSjet_wp, UncertaintySource unc_source,
                                     UncertaintyScale unc_scale)
{
    double weight = 1.;
    for(size_t leg_id = 1; leg_id <= 2; ++leg_id)
        weight *= GetLegIdIsoWeight(eventInfo.GetLeg(leg_id), VSe_wp, VSmu_wp, VSjet_wp, unc_source, unc_scale) *
                  GetCustomTauSF(eventInfo.GetLeg(leg_id), unc_source, unc_scale, eventInfo.GetChannel());
    return weight;
}

double LeptonWeights::GetTriggerWeight(EventInfo& eventInfo, DiscriminatorWP VSjet_wp, UncertaintySource unc_source,
                                       UncertaintyScale unc_scale)
{
    bool data_same_as_central, mc_same_as_central;
    const double eff_data = GetTriggerEfficiency(eventInfo, true, VSjet_wp, unc_source, unc_scale,
                                                 data_same_as_central);
    const double eff_mc = GetTriggerEfficiency(eventInfo, false, VSjet_wp, unc_source, unc_scale, mc_same_as_central);
    if(data_same_as_central && mc_same_as_central)
        return eff_data / eff_mc;

    bool tmp;
    const double eff_data_central = GetTriggerEfficiency(eventInfo, true, VSjet_wp, UncertaintySource::None,
                                                         UncertaintyScale::Central, tmp);
    const double eff_mc_central = GetTriggerEfficiency(eventInfo, false, VSjet_wp, UncertaintySource::None,
                                                       UncertaintyScale::Central, tmp);

    const auto unc_scale_inv = unc_scale == UncertaintyScale::Up ? UncertaintyScale::Down : UncertaintyScale::Up;
    const double eff_data_inv = GetTriggerEfficiency(eventInfo, true, VSjet_wp, unc_source, unc_scale_inv, tmp);
    const double eff_mc_inv = GetTriggerEfficiency(eventInfo, false, VSjet_wp, unc_source, unc_scale_inv, tmp);

    const double err_data = std::max(std::abs(eff_data_central - eff_data), std::abs(eff_data_central - eff_data_inv));
    const double err_mc = std::max(std::abs(eff_mc_central - eff_mc), std::abs(eff_mc_central - eff_mc_inv));

    const PhysicalValue pv_eff_data(eff_data_central, err_data);
    const PhysicalValue pv_eff_mc(eff_mc_central, err_mc);
    const PhysicalValue sf = pv_eff_data / pv_eff_mc;
    return sf.GetValue() + static_cast<int>(unc_scale) * sf.GetStatisticalError();
}

double LeptonWeights::GetTriggerPrescaleWeight(EventInfo& eventInfo) const
{
    //values calculated using brilcalc (https://docs.google.com/spreadsheets/d/1vlca_ap3Ppc62BzjjMcyRNAwUvNhulIU/edit#gid=1424643624)
    //2016
    static const std::vector<std::string> triggerPaths_Prescaled_IsoMu22_2016 = {"HLT_IsoMu22_v", "HLT_IsoMu22_v"};
    static const std::vector<std::string> triggerPaths_Prescaled_IsoMu22_eta2p1_2016 = {
        "HLT_IsoMu22_eta2p1_v", "HLT_IsoTkMu22_eta2p1_v"
    };
    static const std::vector<std::string> triggerPaths_unPrescaled_2016 = {
            "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v" };

    static constexpr double prescaled_weight_IsoMu22_2016 = 0.7952219705;
    static constexpr double prescaled_weight_IsoMu22_eta2p1_2016 = 0.9237823527;

    //2017
    static const std::vector<std::string> triggerPaths_Prescaled_muTau_2017 = { "HLT_IsoMu24_v" };
    static const std::vector<std::string> triggerPaths_unPrescaled_2017 = {
        "HLT_IsoMu27_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" };
    static const std::vector<std::string> triggerPaths_single_unPrescaled_2017 = { "HLT_IsoMu27_v" };

    static constexpr double prescaled_weight_muTau_2017 = 0.9161390096;

    static constexpr double prescaled_weight_VBF_2017 = 0.65308574;
    static constexpr double prescaled_weight_VBF_2018 = 0.990342;

    const Channel channel = eventInfo.GetChannel();
    if(channel == Channel::MuTau || channel == Channel::MuMu) {
        if(period == Period::Run2016) {
            if(channel == Channel::MuTau
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2016))
                return 1.;
            if(eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_IsoMu22_eta2p1_2016))
                return prescaled_weight_IsoMu22_eta2p1_2016;
            return prescaled_weight_IsoMu22_2016;
        } else if(period == Period::Run2017) {
            const auto& paths = channel == Channel::MuTau ? triggerPaths_unPrescaled_2017
                                                          : triggerPaths_single_unPrescaled_2017;
            if(!eventInfo.GetTriggerResults().AnyAcceptAndMatch(paths))
                return prescaled_weight_muTau_2017;
        }
    } else if(channel == Channel::TauTau) {
        if(period == Period::Run2017 && !eventInfo.PassNormalTriggers())
            return prescaled_weight_VBF_2017;
        else if(period == Period::Run2018 && !eventInfo.PassNormalTriggers())
            return prescaled_weight_VBF_2018;
    }
    return 1.;
}

double LeptonWeights::Get(EventInfo& /*eventInfo*/) const
{
    throw exception("Function LeptonWeights::Get no longer suits our needs.");
}

double LeptonWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
}

double LeptonWeights::GetTriggerEfficiency(EventInfo& eventInfo, bool isData, DiscriminatorWP VSjet_wp,
                                           UncertaintySource unc_source, UncertaintyScale unc_scale,
                                           bool& same_as_central)
{
    static const std::map<UncertaintySource, int> tau_dm_unc = {
        { UncertaintySource::TauTriggerUnc, -1 },
        { UncertaintySource::TauTriggerUnc_DM0, 0 },
        { UncertaintySource::TauTriggerUnc_DM1, 1 },
        { UncertaintySource::TauTriggerUnc_DM10, 10 },
        { UncertaintySource::TauTriggerUnc_DM11, 11 },
    };

    static const std::map<UncertaintySource, int> vbf_tau_dm_unc = {
        { UncertaintySource::VBFTauTriggerUnc_DM0, 0 },
        { UncertaintySource::VBFTauTriggerUnc_DM1, 1 },
        { UncertaintySource::VBFTauTriggerUnc_3prong, 10 },
        { UncertaintySource::VBFTauTriggerUnc_3prong, 11 },
    };

    static const std::map<std::pair<Period, bool>, double> ele_min_pt = {
        { { Period::Run2016, false }, 26. },
        { { Period::Run2017, false }, 33. },
        { { Period::Run2017, true }, 25. },
        { { Period::Run2018, false }, 33. },
        { { Period::Run2018, true }, 25. },
    };

    static const std::map<std::pair<Period, bool>, double> muon_min_pt = {
        { { Period::Run2016, false }, 23. },
        { { Period::Run2016, true }, 20. },
        { { Period::Run2017, false }, 25. },
        { { Period::Run2017, true }, 21. },
        { { Period::Run2018, false }, 25. },
        { { Period::Run2018, true }, 21. },
    };

    static const std::map<std::tuple<Period, Channel, bool>, double> tau_min_pt = {
        { { Period::Run2016, Channel::MuTau, false }, 25. },
        { { Period::Run2016, Channel::TauTau, false }, 40. },
        { { Period::Run2017, Channel::ETau, false }, 35. },
        { { Period::Run2017, Channel::MuTau, false }, 32. },
        { { Period::Run2017, Channel::TauTau, false }, 40. },
        { { Period::Run2017, Channel::TauTau, true }, 25. },
        { { Period::Run2018, Channel::ETau, false }, 35. },
        { { Period::Run2018, Channel::MuTau, false }, 32. },
        { { Period::Run2018, Channel::TauTau, false }, 40. },
        { { Period::Run2018, Channel::TauTau, true }, 25. },
    };

    static constexpr double cross_tau_max_eta = 2.1;

    const Channel channel = eventInfo.GetChannel();
    same_as_central = true;

    const auto getEleEff = [&](bool cross) {
        const auto& p4 = eventInfo.GetLeg(1).GetMomentum();
        const auto pt_iter = ele_min_pt.find(std::make_pair(period, cross));
        if(pt_iter == ele_min_pt.end() || p4.pt() <= pt_iter->second) return 0.;
        if(cross && std::abs(eventInfo.GetLeg(2).GetMomentum().eta()) >= cross_tau_max_eta) return 0.;
        const UncertaintyScale scale = unc_source == UncertaintySource::EleTriggerUnc
                                     ? unc_scale : UncertaintyScale::Central;
        if(scale != UncertaintyScale::Central)
            same_as_central = false;
        return cross ? electronSF.GetTriggerEffCross(p4, isData, scale) : electronSF.GetTriggerEff(p4, isData, scale);
    };

    const auto getMuonEff = [&](bool cross) {
        const auto& p4 = eventInfo.GetLeg(1).GetMomentum();
        const auto pt_iter = muon_min_pt.find(std::make_pair(period, cross));
        if(pt_iter == muon_min_pt.end() || p4.pt() <= pt_iter->second) return 0.;
        if(cross && std::abs(eventInfo.GetLeg(2).GetMomentum().eta()) >= cross_tau_max_eta) return 0.;
        const UncertaintyScale scale = unc_source == UncertaintySource::MuonTriggerUnc
                                     ? unc_scale : UncertaintyScale::Central;
        if(scale != UncertaintyScale::Central)
            same_as_central = false;
        return cross ? muonSF.GetTriggerEffCross(p4, isData, scale) : muonSF.GetTriggerEff(p4, isData, scale);
    };

    const auto  match_selection_vbf = [&](bool vbf) {
        //di-tau Trigger
        const auto& leg_1_pt = eventInfo.GetLeg(1).GetMomentum().pt();
        const auto& leg_2_pt = eventInfo.GetLeg(2).GetMomentum().pt();
        if(!vbf && leg_1_pt > 40 && leg_2_pt > 40) return true;

        // VBF tau trigger
        else if(vbf && (period == Period::Run2017 || period == Period::Run2018) && eventInfo.HasVBFjetPair()) {
            const auto& vbf_jet_1_pt = static_cast<float>(eventInfo.GetVBFJet(1).GetMomentum().pt());
            const auto& vbf_jet_2_pt = static_cast<float>(eventInfo.GetVBFJet(2).GetMomentum().pt());
            const auto& vbf_jets_m = static_cast<float>((eventInfo.GetVBFJet(1).GetMomentum() +
                                                         eventInfo.GetVBFJet(2).GetMomentum()).M());
            if(vbf_jets_m > 800 && vbf_jet_1_pt > 140 && vbf_jet_2_pt > 60 && leg_1_pt > 25 && leg_2_pt > 25) return true;
        }
        return false;
    };

    const auto getTauEff = [&](size_t leg_id, bool vbf) {
        const auto& leg = eventInfo.GetLeg(leg_id);
        const float pt = static_cast<float>(leg.GetMomentum().pt());
        const int dm = leg->decayMode();
        const auto pt_iter = tau_min_pt.find(std::make_tuple(period, channel, vbf));
        if(pt_iter == tau_min_pt.end() || pt <= pt_iter->second) return 0.f;
        UncertaintyScale current_scale = UncertaintyScale::Central;
        if(!vbf) {
            const auto& sf_provider = GetTauTriggerSFProvider(channel, VSjet_wp, false);

            auto iter = tau_dm_unc.find(unc_source);
            if(iter != tau_dm_unc.end() && (iter->second < 0 || dm == iter->second))
                current_scale = unc_scale;
            if(current_scale != UncertaintyScale::Central)
                same_as_central = false;
            const int scale = static_cast<int>(current_scale);

            return isData ? sf_provider.getEfficiencyData(pt, dm, scale) : sf_provider.getEfficiencyMC(pt, dm, scale);
        } else {
            const auto& sf_provider_vbf = GetTauTriggerSFProvider(channel, VSjet_wp, true);

            auto iter = vbf_tau_dm_unc.find(unc_source);
            if(iter != vbf_tau_dm_unc.end() && (iter->second < 0 || dm == iter->second))
                current_scale = unc_scale;
            if(current_scale != UncertaintyScale::Central)
                same_as_central = false;
            const int scale = static_cast<int>(current_scale);

            return isData ? sf_provider_vbf.getEfficiencyData(pt, dm, scale) : sf_provider_vbf.getEfficiencyMC(pt, dm,
                                                                                                               scale);
        }
    };

    double efficiency = 0;
    if(channel == Channel::ETau || channel == Channel::MuTau) {
        const double single_eff = channel == Channel::ETau ? getEleEff(false) : getMuonEff(false);
        const double cross_eff_lep = channel == Channel::ETau ? getEleEff(true) : getMuonEff(true);
        const double tau_eff = cross_eff_lep > 0 ? getTauEff(2, false) : 0.;
        const double cross_eff = cross_eff_lep * tau_eff;
        const double cmb_eff = single_eff + cross_eff - std::min(single_eff, cross_eff_lep) * tau_eff;
        efficiency = std::min(cmb_eff, 1.);
    } else if(channel == Channel::MuMu) {
        efficiency = getMuonEff(false);
    } else if(channel == Channel::TauTau) {
        if(match_selection_vbf(false))
           efficiency = getTauEff(1, false) * getTauEff(2, false);
        else if(match_selection_vbf(true))
            efficiency = getTauEff(1, true) * getTauEff(2, true) * static_cast<float>(GetVBFTriggerEfficiency(eventInfo,
                isData, unc_source, unc_scale, same_as_central));
        else efficiency = 0.;
    } else
        throw exception("LeptonWeights::GetTriggerEfficiency: channel is not supported");

    if(!same_as_central && efficiency > 1) efficiency = 1;
    if(!same_as_central && efficiency < 0) efficiency = 0;

    if(efficiency < 0 || efficiency > 1)
        throw exception("LeptonWeights::GetTriggerEfficiency: invalid trigger efficiency = %1% for event %2%.")
                        % eventInfo.GetEventId() % efficiency;
    return efficiency;
}

double LeptonWeights::GetVBFTriggerEfficiency(EventInfo& eventInfo, bool isData, UncertaintySource unc_source,
                                              UncertaintyScale unc_scale, bool& same_as_central)
{
    //SF taken from: https://github.com/camendola/VBFTriggerSFs
    const auto vbf_jet_1_pt = static_cast<float>(eventInfo.GetVBFJet(1).GetMomentum().pt());
    const auto vbf_jet_2_pt = static_cast<float>(eventInfo.GetVBFJet(2).GetMomentum().pt());
    const auto vbf_jets_m = static_cast<float>((eventInfo.GetVBFJet(1).GetMomentum() +
                                                eventInfo.GetVBFJet(2).GetMomentum()).M());
    std::string_view input(tauVBFTriggerInput);
    vbf_trigger_provider = std::make_shared<VBFTriggerSFs>(input);

    UncertaintyScale scale =  UncertaintyScale::Central;
    scale = unc_source == UncertaintySource::VBFTriggerUnc_jets ? unc_scale : UncertaintyScale::Central;
    if(scale != UncertaintyScale::Central) same_as_central = false;
    const int scale_value = static_cast<int>(scale);

    const double eff = isData ? vbf_trigger_provider->getJetsEfficiencyData(vbf_jets_m, vbf_jet_1_pt, vbf_jet_2_pt,
        scale_value) : vbf_trigger_provider->getJetsEfficiencyMC(vbf_jets_m, vbf_jet_1_pt, vbf_jet_2_pt, scale_value);
    return eff;
}

double  LeptonWeights::GetCustomTauSF(const LepCandidate& leg, UncertaintySource unc_source, UncertaintyScale unc_scale,
                                      Channel channel)
{
    // taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DoubleHiggsToBBTauTauWorkingLegacyRun2#Custom_Tau_ID_SF
    static const std::map<int, std::map<double, std::map<UncertaintyScale, double>>> tau_sf_error = {
        { 0, { { 1.078, {{UncertaintyScale::Up, 0.034}, {UncertaintyScale::Down, 0.036} }} } },
        { 1, { { 1.112, {{UncertaintyScale::Up, 0.023}, {UncertaintyScale::Down, 0.023} }} } },
        { 10, { { 0.984, {{UncertaintyScale::Up, 0.063}, {UncertaintyScale::Down, 0.067} }} } },
        { 11, { { 0.759, {{UncertaintyScale::Up, 0.178}, {UncertaintyScale::Down, 0.259} }} } }
    };

    static const std::map<int, UncertaintySource> dm_unc_sources = {
         {0, UncertaintySource::TauCustomSF_DM0},
         {1, UncertaintySource::TauCustomSF_DM1},
         {10, UncertaintySource::TauCustomSF_DM10},
         {11, UncertaintySource::TauCustomSF_DM11}
    };

    if(period == Period::Run2017 && channel == Channel::TauTau && leg->leg_type() == LegType::tau
       && leg->gen_match() == GenLeptonMatch::Tau){
        const auto second_map = tau_sf_error.at(leg->decayMode());
        const UncertaintyScale scale = unc_source == dm_unc_sources.at(leg->decayMode())
                                       ? unc_scale : UncertaintyScale::Central;
        double sf = 1.;
        for(const auto& third_map : second_map){
            double sf_value = third_map.first;
            std::map<UncertaintyScale, double> error_map = third_map.second;

            double error = scale == UncertaintyScale::Central ? 0. :
                           static_cast<int>(scale) * error_map.at(scale);
            sf = sf_value + error;
        }
        return sf;
    }
    return 1.;
}

bool LeptonWeights::ApplyIdUncertaintyScale(int decay_mode, double pt, double eta, GenLeptonMatch gen_lepton_match,
                                            LegType leg_type, UncertaintySource unc_source)
{
    static const std::map<UncertaintySource, std::set<int>> unc_source_map_dm = {
        { UncertaintySource::TauVSjetSF_DM0, {0} },
        { UncertaintySource::TauVSjetSF_DM1, {1} },
        { UncertaintySource::TauVSjetSF_3prong, {10, 11} },
    };

    static const std::map<UncertaintySource, std::pair<double, double>> unc_source_map_pt = {
        { UncertaintySource::TauVSjetSF_pt20to25, {20, 25} },
        { UncertaintySource::TauVSjetSF_pt25to30, {25, 30} },
        { UncertaintySource::TauVSjetSF_pt30to35, {30, 35} },
        { UncertaintySource::TauVSjetSF_pt35to40, {35, 40} },
        { UncertaintySource::TauVSjetSF_ptgt40, {40, std::numeric_limits<double>::max()} },
    };

    static const std::map<UncertaintySource, std::pair<double, double>> unc_source_map_eta_ele = {
        { UncertaintySource::TauVSeSF_barrel, {0, 1.460} },
        { UncertaintySource::TauVSeSF_endcap, {1.558, 3} },
    };

    static const std::map<UncertaintySource, std::pair<double, double>> unc_source_map_eta_muon = {
        { UncertaintySource::TauVSmuSF_etaLt0p4, {0, 0.4} },
        { UncertaintySource::TauVSmuSF_eta0p4to0p8, {0.4, 0.8} },
        { UncertaintySource::TauVSmuSF_eta0p8to1p2, {0.8, 1.2} },
        { UncertaintySource::TauVSmuSF_eta1p2to1p7, {1.2, 1.7} },
        { UncertaintySource::TauVSmuSF_etaGt1p7, {1.7, 3} },
    };

    static const std::map<GenLeptonMatch, UncertaintySource> unc_source_map = {
        { GenLeptonMatch::Electron, UncertaintySource::EleIdIsoUnc },
        { GenLeptonMatch::TauElectron, UncertaintySource::EleIdIsoUnc },
        { GenLeptonMatch::Muon, UncertaintySource::MuonIdIsoUnc },
        { GenLeptonMatch::TauMuon, UncertaintySource::MuonIdIsoUnc },
    };

    if(leg_type == LegType::e || leg_type == LegType::mu) {
        const auto iter = unc_source_map.find(gen_lepton_match);
        return iter != unc_source_map.end() && iter->second == unc_source;
    } else if(leg_type == LegType::tau) {
        if(gen_lepton_match == GenLeptonMatch::Tau) {
            if(is_dm_binned) {
                const auto iter = unc_source_map_dm.find(unc_source);
                return iter != unc_source_map_dm.end() && iter->second.count(decay_mode);
            } else {
                const auto iter = unc_source_map_pt.find(unc_source);
                return iter != unc_source_map_pt.end() && pt >= iter->second.first && pt < iter->second.second;
            }
        } else if(gen_lepton_match == GenLeptonMatch::Electron || gen_lepton_match == GenLeptonMatch::TauElectron ||
                  gen_lepton_match == GenLeptonMatch::Muon ||  gen_lepton_match == GenLeptonMatch::TauMuon) {
            const bool is_ele = gen_lepton_match == GenLeptonMatch::Electron
                                || gen_lepton_match == GenLeptonMatch::TauElectron;
            const auto& eta_map = is_ele ? unc_source_map_eta_ele : unc_source_map_eta_muon;
            const double abs_eta = std::abs(eta);
            const auto iter = eta_map.find(unc_source);
            return iter != eta_map.end() && abs_eta >= iter->second.first && abs_eta < iter->second.second;
        }
        return false;
    } else {
        throw exception("LeptonWeights::ApplyIdUncertaintyScale: unknown leg type");
    }
}

} // namespace mc_corrections
} // namespace analysis
