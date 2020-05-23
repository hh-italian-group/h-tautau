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

bool LeptonScaleFactors::HasCrossTriggers() const { return triggerCross.get() != nullptr; }

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
                             const std::string& tauTriggerInput, Period _period, DiscriminatorWP _tau_iso_wp,  bool _is_dm_binned) :
    electronSF(electron_idIsoInput, electron_SingletriggerInput, electron_CrossTriggerInput),
    muonSF(muon_idIsoInput, muon_SingletriggerInput, muon_CrossTriggerInput), period(_period),
    tau_iso_wp(_tau_iso_wp), is_dm_binned(_is_dm_binned)
{

    tauTriggerWeight_eTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "etau", ToString(tau_iso_wp));
    tauTriggerWeight_muTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "mutau", ToString(tau_iso_wp));
    tauTriggerWeight_tauTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "ditau",
                                                                        ToString(tau_iso_wp));
}

TauIDSFTool& LeptonWeights::GetTauIdProvider(TauIdDiscriminator discr, DiscriminatorWP wp)
{
    static const std::map<Period ,std::string>  period_label = { { Period::Run2016, "2016Legacy" },
                                                     { Period::Run2017, "2017ReReco" },
                                                     { Period::Run2018, "2018ReReco" } };
     auto& provider = tau_sf_providers[discr][wp];
     if(!provider) {
         static const std::string s = "by";
         const auto i = ToString(discr).find(s);
         const std::string discr_name = ToString(discr).erase(i, s.length());
         if(!period_label.count(period))
             throw exception("LeptonWeights::GetTauIdProvider: Period %1% is not supported.") % period;
         provider = std::make_shared<TauIDSFTool>(period_label.at(period), discr_name, ToString(wp), is_dm_binned, false);
     }
     return *provider;
}

double LeptonWeights::GetLegIdIsoWeight(LepCandidate leg, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                                        DiscriminatorWP VSjet_wp,
                                        UncertaintySource unc_source, UncertaintyScale unc_scale)
{
    const UncertaintyScale current_scale = ApplyUncertaintyScale(leg->decayMode(), leg.GetMomentum().pt(),
                                                                 leg.GetMomentum().eta(),leg->gen_match(),
                                                                 leg->leg_type(), unc_source)
                                           ? unc_scale : UncertaintyScale::Central;

    if(leg->leg_type() == LegType::e)
        return electronSF.GetIdIsoSF(leg.GetMomentum()) +
            static_cast<int>(current_scale) * electronSF.GetIdIsoSFError(leg.GetMomentum());

    else if(leg->leg_type() == LegType::mu)
        return  muonSF.GetIdIsoSF(leg.GetMomentum()) +
            static_cast<int>(current_scale) * muonSF.GetIdIsoSFError(leg.GetMomentum());

    else if(leg->leg_type() == LegType::tau){
        double tau_id_weight_vs_mu = 1;
        double tau_id_weight_vs_ele = 1;
        double tau_id_weight = 1;

        const std::string scale = current_scale == UncertaintyScale::Central ? "" : ToString(current_scale);

        if(leg->gen_match() == GenLeptonMatch::Electron || leg->gen_match() == GenLeptonMatch::TauElectron){
            auto tauIdWeightVsEle = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSe, VSe_wp);
            tau_id_weight_vs_ele = tauIdWeightVsEle.getSFvsEta(leg.GetMomentum().eta(),
                                                              static_cast<int>(leg->gen_match()), scale);
        }
        else if(leg->gen_match() == GenLeptonMatch::Muon || leg->gen_match() == GenLeptonMatch::TauMuon ){
            auto tauIdWeightVsMu = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSmu, VSmu_wp);
            tau_id_weight_vs_mu = tauIdWeightVsMu.getSFvsEta(leg.GetMomentum().eta(), static_cast<int>(leg->gen_match()),
                                                             scale);
        }
        else if(leg->gen_match() == GenLeptonMatch::Tau){
            auto tauIdWeight = GetTauIdProvider(TauIdDiscriminator::byDeepTau2017v2p1VSjet, VSjet_wp);
            tau_id_weight = is_dm_binned ? tauIdWeight.getSFvsDM(leg.GetMomentum().pt(), leg->decayMode(),
                                                                        static_cast<int>(leg->gen_match()), scale)
                                                : tauIdWeight.getSFvsPT(leg.GetMomentum().pt(),
                                                                        static_cast<int>(leg->gen_match()), scale);
        }

        return tau_id_weight * tau_id_weight_vs_mu * tau_id_weight_vs_ele;
    }
    else
        throw exception("Leg type not found.");
}

double LeptonWeights::GetIdIsoWeight(EventInfo& eventInfo, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                                     DiscriminatorWP VSjet_wp, UncertaintySource unc_source,
                                     UncertaintyScale unc_scale)
{
    double weight = 1;
    for(size_t leg_id = 1; leg_id <= 2; ++leg_id)
        weight *= GetLegIdIsoWeight(eventInfo.GetLeg(leg_id), VSe_wp, VSmu_wp, VSjet_wp, unc_source, unc_scale);
    return weight;
}

double LeptonWeights::GetTriggerWeight(EventInfo& eventInfo) const
{
    const ntuple::Event& event = *eventInfo;
    const double eff_data = GetTriggerEfficiency(eventInfo, true);
    const double eff_mc = GetTriggerEfficiency(eventInfo, false);
    if(eff_mc == 0) {
        if(event.lep_decayMode.at(0) >= 0 && eff_data != 0)
            throw exception("Undefined trigger SF for event:%1%. eff mc = 0 & eff data= %2%")
                  % EventIdentifier(event) % eff_data;
        return 0;
    }
    return eff_data / eff_mc;
}


double LeptonWeights::Get(EventInfo& /*eventInfo*/) const
{
    throw exception("Function LeptonWeights::Get no longer suits our needs.");
}

double LeptonWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
}

double LeptonWeights::GetTriggerEfficiency(EventInfo& eventInfo, bool isData) const
{
    const Event& event = *eventInfo;
    const Channel channel = static_cast<Channel>(event.channelId);
    double prescaled_weight = 1;

    //values calculated using brilcalc (https://docs.google.com/spreadsheets/d/1vlca_ap3Ppc62BzjjMcyRNAwUvNhulIU/edit#gid=1424643624)
    //2016
    static const std::vector<std::string> triggerPaths_Prescaled_IsoMu22_2016 = {"HLT_IsoMu22_v", "HLT_IsoMu22_v"};
    static const std::vector<std::string> triggerPaths_Prescaled_IsoMu22_eta2p1_2016 = {"HLT_IsoMu22_eta2p1_v", "HLT_IsoTkMu22_eta2p1_v"};
    static constexpr double prescaled_weight_IsoMu22_2016 = 0.7952219705;
    static constexpr double prescaled_weight_IsoMu22_eta2p1_2016 = 0.9237823527;

    //2017
    static const std::vector<std::string> triggerPaths_Prescaled_muTau_2017 = {"HLT_IsoMu24_v"};
    static constexpr double prescaled_weight_muTau_2017 = 0.9161390096;

    if(channel == Channel::ETau) {
        if(electronSF.HasCrossTriggers() && std::abs(eventInfo.GetLeg(2).GetMomentum().eta()) < 2.1){
            const double ele_single_eff = electronSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
            const double ele_cross_eff = electronSF.GetTriggerEffCross(eventInfo.GetLeg(1).GetMomentum(), isData);

            double tau_eff = 1;
            tau_eff = isData ?
                tauTriggerWeight_eTau->getEfficiencyData(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                            eventInfo.GetLeg(2)->decayMode()) :
                tauTriggerWeight_eTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                          eventInfo.GetLeg(2)->decayMode());

            static const std::vector<std::string> triggerPaths_unPrescaled_2017 = {
                "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v" };

            return prescaled_weight * std::min(ele_single_eff * (1 - tau_eff) + ele_cross_eff * tau_eff, 1.);
        }
        else {
            return prescaled_weight * electronSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
        }

    }
    else if(channel == Channel::MuTau) {
        if(muonSF.HasCrossTriggers() && std::abs(eventInfo.GetLeg(2).GetMomentum().eta()) < 2.1){
                const double muon_single_eff = muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
                const double muon_cross_eff = muonSF.GetTriggerEffCross(eventInfo.GetLeg(1).GetMomentum(), isData);
                double tau_eff = 1;
                tau_eff = isData ?
                    tauTriggerWeight_muTau->getEfficiencyData(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                              eventInfo.GetLeg(2)->decayMode()) :
                    tauTriggerWeight_muTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                            eventInfo.GetLeg(2)->decayMode());

                static const std::vector<std::string> triggerPaths_unPrescaled_2016 = {
                        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v" };

                static const std::vector<std::string> triggerPaths_unPrescaled_2017 = {
                    "HLT_IsoMu27_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" };

                if(period == Period::Run2016 &&
                    !(eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2016))
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_IsoMu22_eta2p1_2016))
                    prescaled_weight = prescaled_weight_IsoMu22_eta2p1_2016;

                else if(period == Period::Run2016 &&
                    !(eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2016) ||
                        eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_IsoMu22_eta2p1_2016))
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_IsoMu22_2016))
                    prescaled_weight = prescaled_weight_IsoMu22_2016;

                if(period == Period::Run2017 && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2017)
                        && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau_2017))
                    prescaled_weight = prescaled_weight_muTau_2017;

                return prescaled_weight * std::min(muon_single_eff * (1 - tau_eff) + muon_cross_eff * tau_eff, 1.);
        }
        else{
            static const std::vector<std::string> triggerPaths_unPrescaled_2017 = {
                "HLT_IsoMu27_v"};
            if(period == Period::Run2017
                    && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2017)
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau_2017))
                prescaled_weight = prescaled_weight_muTau_2017;
            return prescaled_weight * muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
        }

    }

    else if(channel == Channel::MuMu){
        //muMu for period Run2016 prescales should also be applied if all single mu triggers are used.
        static const std::vector<std::string> triggerPaths_unPrescaled_2017 = {
            "HLT_IsoMu27_v"};
        if(period == Period::Run2017
                && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled_2017)
                && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau_2017))
            prescaled_weight = prescaled_weight_muTau_2017;
        return  prescaled_weight * muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
    }


    else if(channel == Channel::TauTau){
        double tauSF_1 = 1;
        tauSF_1 = isData ?
            tauTriggerWeight_tauTau->getEfficiencyData(static_cast<float>(eventInfo.GetLeg(1).GetMomentum().pt()),
                                                                          eventInfo.GetLeg(1)->decayMode()) :
            tauTriggerWeight_tauTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(1).GetMomentum().pt()),
                                                                        eventInfo.GetLeg(1)->decayMode());

        double tauSF_2 = 1;
        tauSF_2 = isData ?
            tauTriggerWeight_tauTau->getEfficiencyData(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                          eventInfo.GetLeg(2)->decayMode()) :
            tauTriggerWeight_tauTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                        eventInfo.GetLeg(2)->decayMode());

        return  tauSF_1 * tauSF_2;
    }

    throw exception ("channel not allowed");
}

bool LeptonWeights::ApplyUncertaintyScale(int decay_mode, double pt, double eta, GenLeptonMatch gen_lepton_match,
                                          LegType leg_type, UncertaintySource unc_source)
{
    static const std::map<GenLeptonMatch, std::map<UncertaintySource, std::set<int>>> unc_source_map_dm = {
         {GenLeptonMatch::Tau, { {  UncertaintySource::TauVSjetSF_DM0, {0} } ,
                                 {   UncertaintySource::TauVSjetSF_DM1, {1} } ,
                                 {  UncertaintySource::TauVSjetSF_3prong, {10, 11} } } }
    };

    static const std::map<GenLeptonMatch, std::map<UncertaintySource, std::pair<double, double>>> unc_source_map_pt = {
         {GenLeptonMatch::Tau, { {  UncertaintySource::TauVSjetSF_pt20to25, {20, 25} } ,
                                 {   UncertaintySource::TauVSjetSF_pt25to30,{25, 30} } ,
                                 {  UncertaintySource::TauVSjetSF_pt30to35, {30, 35} },
                                 {  UncertaintySource::TauVSjetSF_pt35to40, {35, 40} },
                                 {  UncertaintySource::TauVSjetSF_ptgt40,   {40, std::numeric_limits<double>::max()} } } }
    };

    static const std::map<GenLeptonMatch, std::map<UncertaintySource, std::pair<double, double>>> unc_source_map_eta = {
         {GenLeptonMatch::Electron, { {  UncertaintySource::TauVSeSF_barrel,        {0, 1.460} } ,
                                      {  UncertaintySource::TauVSeSF_endcap,        {1.558, 3} } } },
         {GenLeptonMatch::Muon,     { {  UncertaintySource::TauVSmuSF_etaLt0p4,     {0, 0.4} } ,
                                      {  UncertaintySource::TauVSmuSF_eta0p4to0p8,  {0.4, 0.8} },
                                      {  UncertaintySource::TauVSmuSF_eta0p8to1p2,  {0.8, 1.2} },
                                      {  UncertaintySource::TauVSmuSF_eta1p2to1p7,  {1.2, 1.7} },
                                      {  UncertaintySource::TauVSmuSF_etaGt1p7,     {1.7, 3} } } },
    };


    static const std::map<GenLeptonMatch, std::map<UncertaintySource, int>> unc_source_map = {
        { GenLeptonMatch::Electron,     { { UncertaintySource::EleIdIsoUnc, -1 } } },
        { GenLeptonMatch::TauElectron,  { { UncertaintySource::EleIdIsoUnc, -1 } } },
        { GenLeptonMatch::Muon,         { { UncertaintySource::MuonIdIsoUnc, -1 } } },
        { GenLeptonMatch::TauMuon,      { { UncertaintySource::MuonIdIsoUnc, -1 } } }
    };

    if(leg_type == LegType::e || leg_type == LegType::mu){
        const auto gen_iter = unc_source_map.find(gen_lepton_match);
        if(gen_iter != unc_source_map.end()) {
            const auto source_iter = gen_iter->second.find(unc_source);
            if(source_iter != gen_iter->second.end()) {
                return source_iter->second < 0;
            }
        }
    }

    if(leg_type == LegType::tau){
        if(gen_lepton_match == GenLeptonMatch::Tau){
            if(is_dm_binned){
                const auto gen_iter = unc_source_map_dm.find(gen_lepton_match);
                if(gen_iter != unc_source_map_dm.end()) {
                    const auto source_iter = gen_iter->second.find(unc_source);
                    if(source_iter != gen_iter->second.end()) {
                        return source_iter->second.count(decay_mode);
                    }
                }
            }
            else{
                const auto gen_iter = unc_source_map_pt.find(gen_lepton_match);
                if(gen_iter != unc_source_map_pt.end()) {
                    const auto source_iter = gen_iter->second.find(unc_source);
                    if(source_iter != gen_iter->second.end()) {
                        auto pt_interval = source_iter->second;
                        return pt >= pt_interval.first  &&  pt < pt_interval.second;
                    }
                }
            }
        }
        else if(gen_lepton_match == GenLeptonMatch::Electron || gen_lepton_match == GenLeptonMatch::TauElectron ||
                gen_lepton_match == GenLeptonMatch::Muon ||  gen_lepton_match == GenLeptonMatch::TauMuon){
            const auto gen_iter = unc_source_map_eta.find(gen_lepton_match);
            if(gen_iter != unc_source_map_eta.end()) {
                const auto source_iter = gen_iter->second.find(unc_source);
                if(source_iter != gen_iter->second.end()) {
                    auto eta_interval = source_iter->second;
                    return std::abs(eta) >= eta_interval.first  &&  std::abs(eta) < eta_interval.second;
                }
            }
        }
    }
    return false;
}


} // namespace mc_corrections
} // namespace analysis
