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
                             const std::string& tauTriggerInput, Period _period, DiscriminatorWP _tau_iso_wp) :
    electronSF(electron_idIsoInput, electron_SingletriggerInput, electron_CrossTriggerInput),
    muonSF(muon_idIsoInput, muon_SingletriggerInput, muon_CrossTriggerInput), period(_period),
    tau_iso_wp(_tau_iso_wp)
{

    tauTriggerWeight_eTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "etau", ToString(tau_iso_wp));
    tauTriggerWeight_muTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "mutau", ToString(tau_iso_wp));
    tauTriggerWeight_tauTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "ditau",
                                                                        ToString(tau_iso_wp));
}

std::shared_ptr<TauIDSFTool>& LeptonWeights::GetTauIdProvider(TauIdDiscriminator discr, DiscriminatorWP wp)
{
    std::map<Period ,std::string>  period_label = { { Period::Run2016, "2016Legacy" },
                                                     { Period::Run2017, "2017ReReco" },
                                                     { Period::Run2018, "2018ReReco" } };
    std::string s = "by";
    std::string::size_type i = ToString(discr).find(s);

    tau_sf_providers[discr][wp] = std::make_shared<TauIDSFTool>(period_label.at(period),
                                ToString(discr).erase(i, s.length()), ToString(wp), false, false);
    return tau_sf_providers[discr][wp];
}

double LeptonWeights::GetLegIdIsoWeight(LepCandidate leg, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                                        DiscriminatorWP VSjet_wp,
                                        UncertaintySource /*unc_source*/, UncertaintyScale /*unc_scale*/)
{
    std::map<LegType, std::pair<TauIdDiscriminator,DiscriminatorWP>> leg_wp_disc = {
         { LegType::e, { TauIdDiscriminator::byDeepTau2017v2p1VSe, VSe_wp } },
         { LegType::mu, { TauIdDiscriminator::byDeepTau2017v2p1VSmu, VSmu_wp } },
         { LegType::tau, { TauIdDiscriminator::byDeepTau2017v2p1VSjet, VSjet_wp } }
     };


    auto tauIdWeight = GetTauIdProvider(leg_wp_disc.at(leg->leg_type()).first, leg_wp_disc.at(leg->leg_type()).second);

    if(leg->leg_type() == LegType::e){
        double tau_id_weight = tauIdWeight->getSFvsEta(leg.GetMomentum().eta(), static_cast<int>(leg->gen_match()));
        return electronSF.GetIdIsoSF(leg.GetMomentum()) * tau_id_weight;
    }
    else if(leg->leg_type() == LegType::mu){
        double tau_id_weight = tauIdWeight->getSFvsEta(leg.GetMomentum().eta(), static_cast<int>(leg->gen_match()));
        return  muonSF.GetIdIsoSF(leg.GetMomentum()) * tau_id_weight;
    }
    else if(leg->leg_type() == LegType::tau){
        double tau_id_weight = tauIdWeight->getSFvsPT(leg.GetMomentum().pt(), static_cast<int>(leg->gen_match()));
        return tau_id_weight;
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

double LeptonWeights::GetTriggerWeight(EventInfo& eventInfo, UncertaintyScale unc_scale) const
{
    const ntuple::Event& event = *eventInfo;
    const double eff_data = GetTriggerEfficiency(eventInfo, true, unc_scale);
    const double eff_mc = GetTriggerEfficiency(eventInfo, false, unc_scale);
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

double LeptonWeights::GetTriggerEfficiency(EventInfo& eventInfo, bool isData, UncertaintyScale unc_scale ) const
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
                                                                            eventInfo.GetLeg(2)->decayMode(),
                                                                            static_cast<int>(unc_scale)) :
                tauTriggerWeight_eTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                          eventInfo.GetLeg(2)->decayMode(),
                                                                          static_cast<int>(unc_scale));

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
                                                              eventInfo.GetLeg(2)->decayMode(),
                                                              static_cast<int>(unc_scale)) :
                    tauTriggerWeight_muTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                            eventInfo.GetLeg(2)->decayMode(),
                                                            static_cast<int>(unc_scale));

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
                                                                          eventInfo.GetLeg(1)->decayMode(),
                                                                          static_cast<int>(unc_scale)) :
            tauTriggerWeight_tauTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(1).GetMomentum().pt()),
                                                                        eventInfo.GetLeg(1)->decayMode(),
                                                                        static_cast<int>(unc_scale));

        double tauSF_2 = 1;
        tauSF_2 = isData ?
            tauTriggerWeight_tauTau->getEfficiencyData(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                          eventInfo.GetLeg(2)->decayMode(),
                                                                          static_cast<int>(unc_scale)) :
            tauTriggerWeight_tauTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                                        eventInfo.GetLeg(2)->decayMode(),
                                                                        static_cast<int>(unc_scale));

        return  tauSF_1 * tauSF_2;
    }

    throw exception ("channel not allowed");
}

} // namespace mc_corrections
} // namespace analysis
