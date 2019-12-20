/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

double TauESUncertainties::GetCorrectionFactor(analysis::Period period, int decayMode, GenLeptonMatch genLeptonMatch,
                                               UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                               TauIdDiscriminator tauVSjetDiscriminator,
                                               TauIdDiscriminator tauVSeDiscriminator, double eta)
{
    if(genLeptonMatch == GenLeptonMatch::Tau) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorTrueTau(period, decayMode, current_scale, pt, tauVSjetDiscriminator);
    }
    else if(genLeptonMatch == GenLeptonMatch::Electron || genLeptonMatch == GenLeptonMatch::TauElectron) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::EleFakingTauES
                                       ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorEleFakingTau(period, current_scale, eta, tauVSeDiscriminator, decayMode);
    }
    else
        return 1.;
}


double TauESUncertainties::GetCorrectionFactorTrueTau(analysis::Period period, int decayMode,
                                                      UncertaintyScale current_scale, double pt,
                                                      TauIdDiscriminator tauVSjetDiscriminator)
{
    //put no shift and 2% of unc for missing DM
    //values taken from: https://indico.cern.ch/event/864131/contributions/3644021/attachments/1946837/3230164/Izaak_TauPOG_TauES_20191118.pdf

    // Corrections for DeepTau
    const static std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor_deep_tau = {
      { analysis::Period::Run2016, { {0, PhysicalValue(-0.1,0.7)},
                                     {1, PhysicalValue(-0.1,0.3)},
                                     {2, PhysicalValue(0.0,2.0)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.0,0.4)},
                                     {11, PhysicalValue(2.6,0.6)}    }},
      { analysis::Period::Run2017, { {0, PhysicalValue(-0.7,0.7)},
                                     {1, PhysicalValue(-1.1,0.3)},
                                     {2, PhysicalValue(0.0,2.0)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.5,0.5)},
                                     {11, PhysicalValue(1.7,0.6)}} },
     { analysis::Period::Run2018,  { {0, PhysicalValue(-1.6,0.8)},
                                     {1, PhysicalValue(0.8,0.3)},
                                     {2, PhysicalValue(0.0,2.0)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(-0.9,0.4)},
                                     {11, PhysicalValue(1.3,1.0)}} }
    };

    // Corrections for MVA
    const static std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor_mva = {
      { analysis::Period::Run2016, { {0, PhysicalValue(-0.6,1.0)},
                                     {1, PhysicalValue(-0.5,0.9)},
                                     {2, PhysicalValue(0.0,2.0)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.0,1.1)},
                                     {11, PhysicalValue(0.0,2.0)}    }},
      { analysis::Period::Run2017, { {0, PhysicalValue(0.7,0.8)},
                                     {1, PhysicalValue(-0.2,0.8)},
                                     {2, PhysicalValue(0.0,2.0)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.1,0.9)},
                                     {11, PhysicalValue(-0.1,1.0)}} },
     { analysis::Period::Run2018, { {0, PhysicalValue(-1.3,1.1)},
                                    {1, PhysicalValue(-0.5,0.9)},
                                    {2, PhysicalValue(0.0,2.0)},
                                    {5, PhysicalValue(0.0,2.0)},
                                    {6, PhysicalValue(0.0,2.0)},
                                    {10, PhysicalValue(-1.2,0.8)},
                                    {11, PhysicalValue(0.0,2.0)}} }
    };

    PhysicalValue tau_correction = PhysicalValue(0.0,0.0);
    if(pt > 400)
        tau_correction = PhysicalValue(0., 3.);
    else {
        const std::map<analysis::Period, std::map<int, PhysicalValue>>* tau_correction_factor = nullptr;

        if(tauVSjetDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSjet)
            tau_correction_factor = &tau_correction_factor_deep_tau;
        else if(tauVSjetDiscriminator == TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017)
            tau_correction_factor = &tau_correction_factor_mva;
        else
            throw analysis::exception("TauIdDiscriminator: '%1%' not allowed.") % tauVSjetDiscriminator;

        if(!tau_correction_factor->count(period))
            throw exception("Period '%1%'not found in tau correction map.") %period;
        if(!tau_correction_factor->at(period).count(decayMode))
            throw exception("Decay mode: '%1%' not found in tau correction map.") %decayMode;

        tau_correction = tau_correction_factor->at(period).at(decayMode);

        //double uncertainty = 1 + ((static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
    }
    auto tau_final_correction = (tau_correction.GetValue() +
                                 static_cast<int>(current_scale) * tau_correction.GetStatisticalError())/100;

    return 1 + tau_final_correction;
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
} // namespace analysis
