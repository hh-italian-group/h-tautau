/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

    double TauESUncertainties::GetCorrectionFactor(analysis::Period period, int decayMode, UncertaintyScale scale,
                                                   double pt, TauIdDiscriminator tauIdDiscriminator)
    {
        //put no shift and 2% of unc for missing DM
        //values taken from: https://indico.cern.ch/event/864131/contributions/3644021/attachments/1946837/3230164/Izaak_TauPOG_TauES_20191118.pdf
        static const std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor_deep_tau = {
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
         { analysis::Period::Run2018, { {0, PhysicalValue(-1.6,0.8)},
                                        {1, PhysicalValue(0.8,0.3)},
                                        {5, PhysicalValue(0.0,2.0)},
                                        {6, PhysicalValue(0.0,2.0)},
                                        {10, PhysicalValue(-0.9,0.4)},
                                        {11, PhysicalValue(1.3,1.0)}} }
        };

        static const std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor_mva = {
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
        // Values taken from: https://indico.cern.ch/event/865792/contributions/3659828/attachments/1954858/3246751/ETauFR-update2Dec.pdf#page=40
        // For eta < 1.448
        static const std::map<DiscriminatorWP, double> deep_tau_vs_e_energy_scale_eta_low = {
          { DiscriminatorWP::VVLoose,  PhysicalValue(1.028,0.003) },
          { DiscriminatorWP::VLoose,  PhysicalValue(1.040,0.003) },
          { DiscriminatorWP::VLoose,  PhysicalValue(1.048,0.005) },
          { DiscriminatorWP::Medium,  PhysicalValue(1.035,0.008) },
          { DiscriminatorWP::VVLoose,  PhysicalValue(1.04,0.03) },
          { DiscriminatorWP::Tight,  PhysicalValue(1.05,0.03) },
          { DiscriminatorWP::VVTight,  PhysicalValue(1.06,0.06) }
        };

        // For eta > 1.558
        static const std::map<DiscriminatorWP, double> deep_tau_vs_e_energy_scale_eta_low = {
          { DiscriminatorWP::VVLoose,  PhysicalValue(0.998,0.003) },
          { DiscriminatorWP::VLoose,  PhysicalValue(0.995,0.005) },
          { DiscriminatorWP::VLoose,  PhysicalValue(1.000,0.015) },
          { DiscriminatorWP::Medium,  PhysicalValue(1.02,0.003) },
          { DiscriminatorWP::VVLoose,  PhysicalValue(1.03,0.06) },
          { DiscriminatorWP::Tight,  PhysicalValue(1.15,0.08) },
          { DiscriminatorWP::VVTight,  PhysicalValue(1.05,0.21) }
        };


        double correction_factor;
        if(pt > 400)
            correction_factor = 1 + static_cast<int>(scale) * 0.03;
        else{
            if(tauIdDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSjet)
                auto tau_correction_factor = tau_correction_factor_deep_tau;
            else if(tauIdDiscriminator == TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017)
                auto tau_correction_factor = tau_correction_factor_mva;
            else
                throw analysis::exception("TauIdDiscriminator: '%1%' not allowed.") % tauIdDiscriminator;

            if(!tau_correction_factor.count(period))
                throw exception("Period not found in tau correction map.");
            if(!tau_correction_factor.at(period).count(decayMode))
                throw exception("Decay mode not found in tau correction map.");
            PhysicalValue tau_correction = tau_correction_factor.at(period).at(decayMode);
            correction_factor = 1 + ((tau_correction.GetValue() + static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
            //double uncertainty = 1 + ((static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
        }
        return correction_factor;
    }

} // namespace analysis
