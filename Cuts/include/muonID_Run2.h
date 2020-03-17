/*! muon ID recommended by the Muon POG.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
    namespace muonID_Run2 {
        constexpr double PFIsoVeryLoose = 0.4; // <
        constexpr double PFIsoLoose = 0.25; // <
        constexpr double PFIsoMedium = 0.20; // <
        constexpr double PFIsoTight = 0.15; // <
        constexpr double PFIsoVeryTight = 0.10; // <
        constexpr double PFIsoVeryVeryTight = 0.05; // <
    }
}
