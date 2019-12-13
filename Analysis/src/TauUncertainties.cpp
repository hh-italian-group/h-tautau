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
                                               TauIdDiscriminator tauVSeDiscriminator, double eta,
                                               DiscriminatorWP tauVSeWP)
{
    if(genLeptonMatch == GenLeptonMatch::Tau) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorTrueTau(period, decayMode, current_scale, pt, tauVSjetDiscriminator);
    }
    else if(genLeptonMatch == GenLeptonMatch::Electron || genLeptonMatch == GenLeptonMatch::TauElectron) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::EleFakingTauES
                                       ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorEleFakingTau(period, current_scale, eta, tauVSeDiscriminator, tauVSeWP);
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
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.0,0.4)},
                                     {11, PhysicalValue(2.6,0.6)}    }},
      { analysis::Period::Run2017, { {0, PhysicalValue(-0.7,0.7)},
                                     {1, PhysicalValue(-1.1,0.3)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.5,0.5)},
                                     {11, PhysicalValue(1.7,0.6)}} },
     { analysis::Period::Run2018,  { {0, PhysicalValue(-1.6,0.8)},
                                     {1, PhysicalValue(0.8,0.3)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(-0.9,0.4)},
                                     {11, PhysicalValue(1.3,1.0)}} }
    };

    // Corrections for MVA
    const static std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor_mva = {
      { analysis::Period::Run2016, { {0, PhysicalValue(-0.6,1.0)},
                                     {1, PhysicalValue(-0.5,0.9)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.0,1.1)},
                                     {11, PhysicalValue(0.0,2.0)}    }},
      { analysis::Period::Run2017, { {0, PhysicalValue(0.7,0.8)},
                                     {1, PhysicalValue(-0.2,0.8)},
                                     {5, PhysicalValue(0.0,2.0)},
                                     {6, PhysicalValue(0.0,2.0)},
                                     {10, PhysicalValue(0.1,0.9)},
                                     {11, PhysicalValue(-0.1,1.0)}} },
     { analysis::Period::Run2018, { {0, PhysicalValue(-1.3,1.1)},
                                    {1, PhysicalValue(-0.5,0.9)},
                                    {5, PhysicalValue(0.0,2.0)},
                                    {6, PhysicalValue(0.0,2.0)},
                                    {10, PhysicalValue(-1.2,0.8)},
                                    {11, PhysicalValue(0.0,2.0)}} }
    };

    double correction_factor = 0;
    PhysicalValue tau_correction = PhysicalValue(0.0,0.0);

    if(pt > 400)
        correction_factor = 1 + static_cast<int>(current_scale) * 0.03;
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
            throw exception("Decay mode not found in tau correction map.");

        tau_correction = tau_correction_factor->at(period).at(decayMode);

        //double uncertainty = 1 + ((static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
    }
    auto tau_final_correction = (tau_correction.GetValue() +
                                 static_cast<int>(current_scale) * tau_correction.GetStatisticalError())/100;

    correction_factor = 1 + tau_final_correction;
    return correction_factor;
}

double TauESUncertainties::GetCorrectionFactorEleFakingTau(analysis::Period period, UncertaintyScale scale, double eta,
                                                           TauIdDiscriminator tauVSeDiscriminator,
                                                           DiscriminatorWP tauVSeWP)
{
    // Values taken from: https://indico.cern.ch/event/865792/contributions/3659828/attachments/1954858/3246751/ETauFR-update2Dec.pdf#page=40
    // For eta < 1.448
    static const std::map<analysis::Period, std::map<DiscriminatorWP, PhysicalValue>> deep_tau_vs_e_energy_scale_eta_low = {
      { analysis::Period::Run2016, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                     { DiscriminatorWP::VVLoose,  PhysicalValue(0.020,0.003)},
                                     { DiscriminatorWP::VLoose,  PhysicalValue(0.030,0.003) },
                                     { DiscriminatorWP::Loose,  PhysicalValue(0.033,0.005) },
                                     { DiscriminatorWP::Medium,  PhysicalValue(0.020,0.015) },
                                     { DiscriminatorWP::Tight,  PhysicalValue(0.06,0.04) },
                                     { DiscriminatorWP::VTight,  PhysicalValue(0.10,0.06) },
                                     { DiscriminatorWP::VVTight,  PhysicalValue(0.03,0.10)}} },
      { analysis::Period::Run2017, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                     { DiscriminatorWP::VVLoose,  PhysicalValue(0.025,0.003)},
                                     { DiscriminatorWP::VLoose,  PhysicalValue(0.035,0.003) },
                                     { DiscriminatorWP::Loose,  PhysicalValue(0.038,0.008) },
                                     { DiscriminatorWP::Medium,  PhysicalValue(0.015,0.018) },
                                     { DiscriminatorWP::Tight,  PhysicalValue(0.04,0.05) },
                                     { DiscriminatorWP::VTight,  PhysicalValue(0.04,0.09) },
                                     { DiscriminatorWP::VVTight,  PhysicalValue(0.07,0.19)}} },
      { analysis::Period::Run2018, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                     { DiscriminatorWP::VVLoose,  PhysicalValue(0.028,0.003)},
                                     { DiscriminatorWP::VLoose,  PhysicalValue(0.040,0.003) },
                                     { DiscriminatorWP::Loose,  PhysicalValue(0.048,0.005) },
                                     { DiscriminatorWP::Medium,  PhysicalValue(0.035,0.008) },
                                     { DiscriminatorWP::Tight,  PhysicalValue(0.04,0.03) },
                                     { DiscriminatorWP::VTight,  PhysicalValue(0.05,0.03) },
                                     { DiscriminatorWP::VVTight,  PhysicalValue(0.06,0.06)}} }
    };

    // For eta > 1.558
    static const std::map<analysis::Period, std::map<DiscriminatorWP, PhysicalValue>> deep_tau_vs_e_energy_scale_eta_high = {
        { analysis::Period::Run2016, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                       { DiscriminatorWP::VVLoose,  PhysicalValue(0.010,0.003)},
                                       { DiscriminatorWP::VLoose,  PhysicalValue(0.023,0.008) },
                                       { DiscriminatorWP::Loose,  PhysicalValue(0.043,0.018) },
                                       { DiscriminatorWP::Medium,  PhysicalValue(0.07,0.04) },
                                       { DiscriminatorWP::Tight,  PhysicalValue(0.10,0.018) },
                                       { DiscriminatorWP::VTight,  PhysicalValue(0.12,0.07) },
                                       { DiscriminatorWP::VVTight,  PhysicalValue(0.11,0.10)}} },
        { analysis::Period::Run2017, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                       { DiscriminatorWP::VVLoose,  PhysicalValue(0.008,0.003)},
                                       { DiscriminatorWP::VLoose,  PhysicalValue(0.018,0.008) },
                                       { DiscriminatorWP::Loose,  PhysicalValue(0.03,0.02) },
                                       { DiscriminatorWP::Medium,  PhysicalValue(0.04,0.04) },
                                       { DiscriminatorWP::Tight,  PhysicalValue(0.11,0.06) },
                                       { DiscriminatorWP::VTight,  PhysicalValue(0.12,0.14) },
                                       { DiscriminatorWP::VVTight,  PhysicalValue(0.08,0.18)}} },
        { analysis::Period::Run2018, { { DiscriminatorWP::VVVLoose,  PhysicalValue(0.0,0.0)},
                                       { DiscriminatorWP::VVLoose,  PhysicalValue(0.998,0.003)},
                                       { DiscriminatorWP::VLoose,  PhysicalValue(0.995,0.005) },
                                       { DiscriminatorWP::Loose,  PhysicalValue(0.0,0.015) },
                                       { DiscriminatorWP::Medium,  PhysicalValue(0.02,0.03) },
                                       { DiscriminatorWP::Tight,  PhysicalValue(0.03,0.06) },
                                       { DiscriminatorWP::VTight,  PhysicalValue(0.15,0.08) },
                                       { DiscriminatorWP::VVTight,  PhysicalValue(0.05,0.21)}} }
      };

    if(!deep_tau_vs_e_energy_scale_eta_high.count(period) || !deep_tau_vs_e_energy_scale_eta_low.count(period) )
        throw exception("Period '%1%'not found in tau vs ele map.") %period;

    PhysicalValue e_fake_rate_correction = PhysicalValue(0.0,0.0);

    if(tauVSeDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSe){

        if(std::abs(eta) < 1.448){
            if(deep_tau_vs_e_energy_scale_eta_low.at(period).count(tauVSeWP))
                e_fake_rate_correction = deep_tau_vs_e_energy_scale_eta_low.at(period).at(tauVSeWP);
            else
                throw exception ("tau vs e wp (low eta map): '%1%' now allowed") %tauVSeWP ;
        }
        else if (std::abs(eta) > 1.558){
            if(deep_tau_vs_e_energy_scale_eta_high.at(period).count(tauVSeWP))
                e_fake_rate_correction = deep_tau_vs_e_energy_scale_eta_high.at(period).at(tauVSeWP);
            else
                throw exception ("tau vs e wp (high eta map): '%1%' now allowed") %tauVSeWP ;
        }

    }

    auto e_fake_rate_final_correction = (e_fake_rate_correction.GetValue() +
                                         static_cast<int>(scale) * e_fake_rate_correction.GetStatisticalError());

    return 1 + e_fake_rate_final_correction;
}
} // namespace analysis
