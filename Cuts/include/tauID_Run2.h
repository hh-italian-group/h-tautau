/*! tau ID recommended by the TAU POG.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
    namespace tauID_Run2 {
        constexpr double pt = 20; // GeV >
        constexpr double eta = 2.3; // <
        constexpr double eta_diTauTrigger = 2.1; // <
        constexpr double dz = 0.2; // <
    }
}
