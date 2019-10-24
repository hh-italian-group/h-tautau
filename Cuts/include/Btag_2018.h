/*! b-jet selection recommended by b-tag POG.
If not specified otherwise, all definitions are taken from the b-tag POG 94X recommendation TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace btag_2018 { 
    constexpr double deepCSVv2L = 0.1241; // >
    constexpr double deepCSVv2M = 0.4184; // >
    constexpr double deepCSVv2T = 0.7527; // >

    constexpr double deepFlavourL = 0.0494; // >
    constexpr double deepFlavourM = 0.2770; // >
    constexpr double deepFlavourT = 0.7264; // >

    constexpr double pt = 20; // >
    constexpr double eta = 2.4; // <
}
} // namespace cuts
