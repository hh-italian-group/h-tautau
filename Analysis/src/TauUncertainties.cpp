/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

double TauESUncertainties::GetCorrectionFactor(analysis::Period period, int decayMode, GenLeptonMatch genLeptonMatch,
                                               UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                               TauIdDiscriminator tauVSeDiscriminator, double eta,
                                               std::string file_tes_low_pt, std::string file_tes_high_pt)
{
    if(genLeptonMatch == GenLeptonMatch::Tau) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorTrueTau(pt, decayMode, file_tes_low_pt, file_tes_high_pt, current_scale,
                                                      genLeptonMatch, unc_source);

    }
    else if(genLeptonMatch == GenLeptonMatch::Electron || genLeptonMatch == GenLeptonMatch::TauElectron) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::EleFakingTauES
                                       ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorEleFakingTau(period, current_scale, eta, tauVSeDiscriminator, decayMode);
    }
    else if(genLeptonMatch == GenLeptonMatch::Muon || genLeptonMatch == GenLeptonMatch::TauMuon)
        return GetCorrectionFactorMuonFakingTau(period, decayMode);
    else
        return 1.;
}


double TauESUncertainties::GetCorrectionFactorTrueTau(double pt, int decayMode,
                                                      std::string file_low_pt, std::string file_high_pt,
                                                      UncertaintyScale scale, GenLeptonMatch genLeptonMatch,
                                                      UncertaintySource unc_source)
{
    auto low_pt = root_ext::OpenRootFile(file_low_pt);
    auto hist_low_pt = std::shared_ptr<TH1F>(root_ext::ReadObject<TH1F>(*low_pt, "tes"));

    auto high_pt = root_ext::OpenRootFile(file_high_pt);
    auto hist_high_pt = std::shared_ptr<TH1F>(root_ext::ReadObject<TH1F>(*high_pt, "tes"));

    std::vector<int> dms = {0, 1, 10, 11};
    if(genLeptonMatch == GenLeptonMatch::Tau && (std::find(dms.begin(), dms.end(), decayMode) != dms.end())){
        if(pt < 100){
            Int_t bin = hist_low_pt->GetXaxis()->FindBin(decayMode);
            double tes  = hist_low_pt->GetBinContent(bin);
            double tes_error = hist_low_pt->GetBinError(bin);

            double tau_final_correction = tes + static_cast<int>(scale) * tes_error;
            return tau_final_correction;
        }
        else {
            Int_t bin = hist_high_pt->GetXaxis()->FindBin(decayMode);
            double tes  = hist_high_pt->GetBinContent(bin);
            double tes_error = hist_high_pt->GetBinError(bin);

            UncertaintyScale current_scale = unc_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
            double tau_final_correction = tes + static_cast<int>(current_scale) * tes_error;
            return tau_final_correction;
        }
    }
    return 1.;
}

double TauESUncertainties::GetCorrectionFactorEleFakingTau(analysis::Period period, UncertaintyScale scale, double eta,
                                                           TauIdDiscriminator tauVSeDiscriminator, int decayMode)
{
    // Values taken from: https://indico.cern.ch/event/868279/contributions/3665970/attachments/1959265/3266335/FES_9Dec_explained.pdf#page=29
    // For eta < 1.448
    static const std::map<analysis::Period, std::map<int, analysis::StVariable>> deep_tau_vs_e_energy_scale_barrel = {
      { analysis::Period::Run2016, { {0,  StVariable(0.679, 0.806, 0.982)},
                                     {1,  StVariable(3.389, 1.168, 2.475)},
                                     {2,  StVariable(3.389, 1.168, 2.475)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }},
      { analysis::Period::Run2017, { {0,  StVariable(0.911, 1.343, 0.882)},
                                     {1,  StVariable(1.154, 2.162, 0.973)},
                                     {2,  StVariable(1.154, 2.162, 0.973)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }},
      { analysis::Period::Run2018, { {0,  StVariable(1.362, 0.904, 0.474)},
                                     {1,  StVariable(1.945, 1.226, 1.598)},
                                     {2,  StVariable(1.945, 1.226, 1.598)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }}
    };

    // For eta > 1.558
    static const std::map<analysis::Period, std::map<int, analysis::StVariable>> deep_tau_vs_e_energy_scale_endcap = {
      { analysis::Period::Run2016, { {0,  StVariable(-3.5, 1.808, 1.102)},
                                     {1,  StVariable(5.0, 5.694, 6.57)},
                                     {2,  StVariable(5.0, 5.694, 6.57)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }},
      { analysis::Period::Run2017, { {0,  StVariable(-2.604, 2.249, 1.43)},
                                     {1,  StVariable(1.5, 4.969, 6.461)},
                                     {2,  StVariable(1.5, 4.969, 6.461)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }},
      { analysis::Period::Run2018, { {0,  StVariable(-3.097, 3.404, 1.25)},
                                     {1,  StVariable(-1.85, 3.772, 5.742)},
                                     {2,  StVariable(-1.85, 3.772, 5.742)},
                                     {5,  StVariable(0.0, 0.0, 0.0)},
                                     {6,  StVariable(0.0, 0.0, 0.0)},
                                     {10, StVariable(0.0, 0.0, 0.0)},
                                     {11, StVariable(0.0, 0.0, 0.0)}    }}
    };

    if(!deep_tau_vs_e_energy_scale_endcap.count(period) || !deep_tau_vs_e_energy_scale_barrel.count(period) )
        throw exception("Period '%1%'not found in tau vs ele map.") %period;

    StVariable e_fake_rate_correction = StVariable(0.0, 0.0, 0.0);

    if(tauVSeDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSe){

        if(std::abs(eta) < 1.460){
            if(deep_tau_vs_e_energy_scale_barrel.at(period).count(decayMode))
                e_fake_rate_correction = deep_tau_vs_e_energy_scale_barrel.at(period).at(decayMode);
            else
                throw exception ("Decay mode for the fake e->tau (barrel): '%1%' now allowed") %decayMode ;
        }
        else if (std::abs(eta) > 1.558){
            if(deep_tau_vs_e_energy_scale_endcap.at(period).count(decayMode))
                e_fake_rate_correction = deep_tau_vs_e_energy_scale_endcap.at(period).at(decayMode);
            else
                throw exception ("Decay mode for the fake e->tau (endcap): '%1%' now allowed") %decayMode ;
        }

    }
    auto error = static_cast<int>(scale) > 0 ? e_fake_rate_correction.error_up : e_fake_rate_correction.error_low;
    auto e_fake_rate_final_correction = (e_fake_rate_correction.value +
                                         static_cast<int>(scale) * error) / 100;

    return 1 + e_fake_rate_final_correction;
}

double TauESUncertainties::GetCorrectionFactorMuonFakingTau(analysis::Period period, int decayMode)
{

    //values taken from: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorkingLegacyRun2#mu_tau_ES
    static const std::map<analysis::Period, std::map<int, float>> mu_correction_factor = {
      { analysis::Period::Run2016, { {0, 0},
                                     {1,  -0.5},
                                     {2,  -0.5},
                                     {5,  0},
                                     {6,  0},
                                     {10, 0},
                                     {11, 0 }}},
      { analysis::Period::Run2017, { {0, -0.2},
                                     {1,  -0.8},
                                     {2,  -0.8},
                                     {5,  0},
                                     {6,  0},
                                     {10, 0},
                                     {11, 0 }}},
      { analysis::Period::Run2018, { {0, -0.2},
                                     {1,  -1.0},
                                     {2,  -1.0},
                                     {5,  0},
                                     {6,  0},
                                     {10, 0},
                                     {11, 0 }}}
    };

    if(!mu_correction_factor.count(period))
        throw exception("Period '%1%'not found in tau correction map.") %period;
    if(!mu_correction_factor.at(period).count(decayMode))
        throw exception("Decay mode: '%1%' not found in tau correction map.") %decayMode;

    const float mu_correction = mu_correction_factor.at(period).at(decayMode);

    return 1 + mu_correction;

}
} // namespace analysis
