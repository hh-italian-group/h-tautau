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
                             const std::string& tauTriggerInput, const std::string& tauTriggerInputOld, Period period,
                             DiscriminatorWP _tau_iso_wp) :
    electronSF(electron_idIsoInput, electron_SingletriggerInput, electron_CrossTriggerInput),
    muonSF(muon_idIsoInput, muon_SingletriggerInput, muon_CrossTriggerInput),
    tau_iso_wp(_tau_iso_wp)
{
    if(period == Period::Run2016){
        tauTriggerWeight =  std::make_shared<TauTriggerWeight2016>(tauTriggerInput);
        tauIdWeight = std::make_shared<TauIdWeight2016>();
    }
    else if(period == Period::Run2017){
        tauTriggerWeight =  std::make_shared<TauTriggerWeight2017>(tauTriggerInput, tauTriggerInputOld,  boost::algorithm::to_lower_copy(ToString(_tau_iso_wp)));
        tauIdWeight = std::make_shared<TauIdWeight2017>();
    }
    else
        throw exception("Period %1% is not supported.") % period;
}

double LeptonWeights::GetIdIsoWeight(const Event& event) const
{
    const Channel channel = static_cast<Channel>(event.channelId);
    if(channel == Channel::ETau) {
        return electronSF.GetIdIsoSF(event.p4_1) * tauIdWeight->GetIdIsoSF(event.p4_2,
            static_cast<GenMatch>(event.gen_match_2), event.decayMode_2, DiscriminatorWP::Tight,
            DiscriminatorWP::Loose, tau_iso_wp);
    }
    else if(channel == Channel::MuTau) {
        return muonSF.GetIdIsoSF(event.p4_1) * tauIdWeight->GetIdIsoSF(event.p4_2,
            static_cast<GenMatch>(event.gen_match_2), event.decayMode_2,  DiscriminatorWP::VLoose,
            DiscriminatorWP::Tight, tau_iso_wp);
    }
    else if(channel == Channel::TauTau) {
        return tauIdWeight->GetIdIsoSF(event.p4_1,
            static_cast<GenMatch>(event.gen_match_1), event.decayMode_1, DiscriminatorWP::VLoose,
            DiscriminatorWP::Loose, tau_iso_wp) * tauIdWeight->GetIdIsoSF(event.p4_2,
            static_cast<GenMatch>(event.gen_match_2), event.decayMode_2,DiscriminatorWP::VLoose,
            DiscriminatorWP::Loose, tau_iso_wp);
    }
    else if(channel == Channel::MuMu)
        return muonSF.GetIdIsoSF(event.p4_1) * muonSF.GetIdIsoSF(event.p4_2);
    else
        throw exception ("channel not allowed");
}

double LeptonWeights::GetTriggerWeight(const Event& event) const
{
    const double eff_data = GetTriggerEfficiency(event, true);
    const double eff_mc = GetTriggerEfficiency(event, false);
    if(eff_mc == 0) {
        if(eff_data != 0)
            throw exception("Undefined trigger SF for %1%") % EventIdentifier(event);
        return 0;
    }
    return eff_data / eff_mc;
}


double LeptonWeights::Get(const Event& event) const { return GetIdIsoWeight(event) * GetTriggerWeight(event); }

double LeptonWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
}

double LeptonWeights::GetTriggerEfficiency(const Event& event, bool isData) const
{
    const Channel channel = static_cast<Channel>(event.channelId);
    if(channel == Channel::ETau) {
        if(electronSF.HasCrossTriggers() && std::abs(event.p4_2.eta()) < 2.1){
            const double ele_single_eff = electronSF.GetTriggerEff(event.p4_1, isData);
            const double ele_cross_eff = electronSF.GetTriggerEffCross(event.p4_1, isData);
            const double tau_eff = tauTriggerWeight->GetEfficiency(channel, LorentzVectorM(event.p4_2),
                                                                   static_cast<GenMatch>(event.gen_match_2),
                                                                   event.decayMode_2, tau_iso_wp, isData);
            return std::min(ele_single_eff * (1 - tau_eff) + ele_cross_eff * tau_eff, 1.);
        }
        else
            return electronSF.GetTriggerEff(event.p4_1, isData);
    }
    else if(channel == Channel::MuTau) {
        if(muonSF.HasCrossTriggers() && std::abs(event.p4_2.eta()) < 2.1){
                const double muon_single_eff = muonSF.GetTriggerEff(event.p4_1, isData);
                const double muon_cross_eff = muonSF.GetTriggerEffCross(event.p4_1, isData);
                const double tau_eff = tauTriggerWeight->GetEfficiency(channel, LorentzVectorM(event.p4_2),
                                                                       static_cast<GenMatch>(event.gen_match_2),
                                                                       event.decayMode_2, tau_iso_wp, isData );
                return std::min(muon_single_eff * (1 - tau_eff) + muon_cross_eff * tau_eff, 1.);
        }
        else
            return muonSF.GetTriggerEff(event.p4_1, isData);
    }

    else if(channel == Channel::MuMu)
        return  muonSF.GetTriggerEff(event.p4_1, isData); // * muonSF.GetTriggerEff(event.p4_2, isData);

    else if(channel == Channel::TauTau){
        double tauSF_1 = tauTriggerWeight->GetEfficiency(channel, LorentzVectorM(event.p4_1), static_cast<GenMatch>(event.gen_match_1),
                                               event.decayMode_1, tau_iso_wp, isData);
        double tauSF_2 = tauTriggerWeight->GetEfficiency(channel, LorentzVectorM(event.p4_2), static_cast<GenMatch>(event.gen_match_2),
                                                   event.decayMode_2, tau_iso_wp, isData);
        return  tauSF_1 * tauSF_2;
    }

    throw exception ("channel not allowed");
}

} // namespace mc_corrections
} // namespace analysis
