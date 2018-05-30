/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace analysis {
namespace uncertainties {
    namespace tau_2016 {
        constexpr double energyUncertainty = 0.012;
        constexpr double sf_1prong = 1 - 0.005;
        constexpr double sf_1prongPi0 = 1 + 0.011;
        constexpr double sf_3prong = 1 + 0.06;
    }

    namespace tau_2017 {
        constexpr double energyUncertainty = 0.03;
        constexpr double sf_1prong = 0.97;
        constexpr double sf_1prongPi0 = 0.98;
        constexpr double sf_3prong = 0.99;
    }

} // namespace uncertainties
} // namespace analysis
