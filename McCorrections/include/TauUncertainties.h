/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {
namespace uncertainties {
    namespace tau_2016 {
        constexpr double energyUncertainty = 0.012;
        constexpr double sf_1prong = 0.994;
        constexpr double sf_1prongPi0 = 0.995;
        constexpr double sf_3prong = 1;
    }

    namespace tau_2017 {
        constexpr double energyUncertainty = 0.03;
        constexpr double sf_1prong = 1.007;
        constexpr double sf_1prongPi0 = 0.998;
        constexpr double sf_3prong = 1.001;
        constexpr double sf_3prongPi0 = 0.999;
    }

    namespace tau_2018 {
        constexpr double energyUncertainty = 0.03;
        constexpr double sf_1prong = 0.987;
        constexpr double sf_1prongPi0 = 0.995;
        constexpr double sf_3prong = 0.988;
    }

} // namespace uncertainties
} // namespace analysis
