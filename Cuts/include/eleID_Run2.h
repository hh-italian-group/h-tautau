/*! electron ID recommended by the e-gamma POG.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
    namespace electronID_Run2 {
        constexpr double high_pt_category = 10; // GeV >

        constexpr char mvaEleID_noIso_Loose[] = "mvaEleID-Fall17-noIso-V2-wpLoose";
        constexpr char mvaEleID_noIso_Medium[] = "mvaEleID-Fall17-noIso-V2-wp90";
        constexpr char mvaEleID_noIso_Tight[] = "mvaEleID-Fall17-noIso-V2-wp80";

        constexpr char mvaEleID_iso_Loose[] = "mvaEleID-Fall17-iso-V2-wpLoose";
        constexpr char mvaEleID_iso_Medium[] = "mvaEleID-Fall17-iso-V2-wp90";
        constexpr char mvaEleID_iso_Tight[] = "mvaEleID-Fall17-iso-V2-wp80";
    }
}
