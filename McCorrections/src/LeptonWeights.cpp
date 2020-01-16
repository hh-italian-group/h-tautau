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
    if(period == Period::Run2016)
        tauIdWeight = std::make_shared<TauIDSFTool>("2016Legacy","DeepTau2017v2p1VSjet", ToString(tau_iso_wp), false);
    else if(period == Period::Run2017)
        tauIdWeight = std::make_shared<TauIDSFTool>("2017ReReco","DeepTau2017v2p1VSjet", ToString(tau_iso_wp), false);
    else if(period == Period::Run2018)
        tauIdWeight = std::make_shared<TauIDSFTool>("2018ReReco","DeepTau2017v2p1VSjet", ToString(tau_iso_wp), false);
    else
        throw exception("Period %1% is not supported.") % period;

    tauTriggerWeight_eTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "etau", ToString(tau_iso_wp));
    tauTriggerWeight_muTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "mutau", ToString(tau_iso_wp));
    tauTriggerWeight_tauTau = std::make_shared<tau_trigger::SFProvider>(tauTriggerInput, "ditau",
                                                                        ToString(tau_iso_wp));
}

double LeptonWeights::GetIdIsoWeight(EventInfoBase& eventInfo) const
{
    const ntuple::Event& event = *eventInfo;
    const Channel channel = static_cast<Channel>(event.channelId);
    if(channel == Channel::ETau) {
        double tau_weight = tauIdWeight->getSFvsPT(eventInfo.GetLeg(2).GetMomentum().pt(),
                                                   static_cast<int>(eventInfo.GetLeg(2)->gen_match()));
        return electronSF.GetIdIsoSF(eventInfo.GetLeg(1).GetMomentum()) * tau_weight;
    }
    else if(channel == Channel::MuTau) {
        double tau_weight = tauIdWeight->getSFvsPT(eventInfo.GetLeg(2).GetMomentum().pt(),
                                                   static_cast<int>(eventInfo.GetLeg(2)->gen_match()));
        return muonSF.GetIdIsoSF(eventInfo.GetLeg(1).GetMomentum()) * tau_weight;
    }
    else if(channel == Channel::TauTau) {
        double tau_weight = tauIdWeight->getSFvsPT(eventInfo.GetLeg(1).GetMomentum().pt(),
                                                   static_cast<int>(eventInfo.GetLeg(1)->gen_match()))
                          * tauIdWeight->getSFvsPT(eventInfo.GetLeg(2).GetMomentum().pt(),
                                                   static_cast<int>(eventInfo.GetLeg(2)->gen_match()));
        return tau_weight;
    }
    else if(channel == Channel::MuMu)
        return muonSF.GetIdIsoSF(eventInfo.GetLeg(1).GetMomentum()) * muonSF.GetIdIsoSF(eventInfo.GetLeg(2).GetMomentum());
    else
        throw exception ("channel not allowed");
}

double LeptonWeights::GetTriggerWeight(EventInfoBase& eventInfo) const
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


double LeptonWeights::Get(EventInfoBase& eventInfo) const
{
    return GetIdIsoWeight(eventInfo) * GetTriggerWeight(eventInfo);
}

double LeptonWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
}

double LeptonWeights::GetTriggerEfficiency(EventInfoBase& eventInfo, bool isData) const
{
    const Event& event = *eventInfo;
    const Channel channel = static_cast<Channel>(event.channelId);
    double prescaled_weight = 1;
    static const std::vector<std::string> triggerPaths_Prescaled_eTau = {"HLT_Ele32_WPTight_Gsf_v"};
    static const std::vector<std::string> triggerPaths_Prescaled_muTau = {"HLT_IsoMu24_v"};
    //Values to be review
    static constexpr double prescaled_weight_eTau = 0.9534869;
    static constexpr double prescaled_weight_muTau = 0.947434;
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

            static const std::vector<std::string> triggerPaths_unPrescaled = {
                "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v" };
            if(period == Period::Run2017
                    && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled)
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_eTau))
                prescaled_weight = prescaled_weight_eTau;
            return prescaled_weight * std::min(ele_single_eff * (1 - tau_eff) + ele_cross_eff * tau_eff, 1.);
        }
        else {
            static const std::vector<std::string> triggerPaths_unPrescaled = {
                "HLT_Ele35_WPTight_Gsf_v"};
            if(period == Period::Run2017
                    && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled)
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_eTau))
                prescaled_weight = prescaled_weight_eTau;
            return electronSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
        }

    }
    else if(channel == Channel::MuTau) {
        if(muonSF.HasCrossTriggers() && std::abs(eventInfo.GetLeg(2).GetMomentum().eta()) < 2.1){
                const double muon_single_eff = muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
                const double muon_cross_eff = muonSF.GetTriggerEffCross(eventInfo.GetLeg(1).GetMomentum(), isData);
                double tau_eff = 1;
                tau_eff = isData ?
                    tauTriggerWeight_muTau->getEfficiencyData(
                            static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                            eventInfo.GetLeg(2)->decayMode()) :
                    tauTriggerWeight_muTau->getEfficiencyMC(static_cast<float>(eventInfo.GetLeg(2).GetMomentum().pt()),
                                                            eventInfo.GetLeg(2)->decayMode());

                static const std::vector<std::string> triggerPaths_unPrescaled = {
                    "HLT_IsoMu27_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" };
                if(period == Period::Run2017
                        && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled)
                        && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau))
                    prescaled_weight = prescaled_weight_muTau;
                return prescaled_weight * std::min(muon_single_eff * (1 - tau_eff) + muon_cross_eff * tau_eff, 1.);
        }
        else{
            static const std::vector<std::string> triggerPaths_unPrescaled = {
                "HLT_IsoMu27_v"};
            if(period == Period::Run2017
                    && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled)
                    && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau))
                prescaled_weight = prescaled_weight_muTau;
            return prescaled_weight * muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
        }

    }

    else if(channel == Channel::MuMu){
        static const std::vector<std::string> triggerPaths_unPrescaled = {
            "HLT_IsoMu27_v"};
        if(period == Period::Run2017
                && !eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_unPrescaled)
                && eventInfo.GetTriggerResults().AnyAcceptAndMatch(triggerPaths_Prescaled_muTau))
            prescaled_weight = prescaled_weight_muTau;
        return  prescaled_weight * muonSF.GetTriggerEff(eventInfo.GetLeg(1).GetMomentum(), isData);
                // * muonSF.GetTriggerEff(event.p4_2, isData);
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

} // namespace mc_corrections
} // namespace analysis
