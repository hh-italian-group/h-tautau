/*! b-jet selection recommended by b-tag POG.
If not specified otherwise, all definitions are taken from the b-tag POG 80X recommendation TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace btag_2016 {
    constexpr double CSVv2L = 0.5426; // >
    constexpr double CSVv2M = 0.8484; // >
    constexpr double CSVv2T = 0.9535; // >

    constexpr double DeepCSVL = 0.2219; // >
    constexpr double DeepCSVM = 0.6324; // >
    constexpr double DeepCSVT = 0.8958; // >

    constexpr double pt = 20; // >
    constexpr double eta = 2.4; // <
}
} // namespace cuts
