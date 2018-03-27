/*! b-jet selection recommended by b-tag POG.
If not specified otherwise, all definitions are taken from the b-tag POG 94X recommendation TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace btag_2017 {
    constexpr double CSVv2L = 0.5803; // >
    constexpr double CSVv2M = 0.8838; // >
    constexpr double CSVv2T = 0.9693; // >
    
    constexpr double deepCSVv2L = 0.1522; // >
    constexpr double deepCSVv2M = 0.4941; // >
    constexpr double deepCSVv2T = 0.8001; // >

    constexpr double pt = 20; // >
    constexpr double eta = 2.4; // <
}
} // namespace cuts
