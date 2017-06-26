/*! Definition of weighting modes for H and HH analyses.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <set>
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {
namespace mc_corrections {

ENUM_OSTREAM_OPERATORS()
ENUM_ISTREAM_OPERATORS()

enum class WeightType {
    PileUp = 0, LeptonTrigIdIso = 1, BTag = 2, DY = 3, TTbar = 4, Wjets = 5, BSM_to_SM = 6, TopPt = 7
};
ENUM_NAMES(WeightType) = {
    { WeightType::PileUp, "PileUp" },
    { WeightType::LeptonTrigIdIso, "LeptonTrigIdIso" },
    { WeightType::BTag, "BTag" },
    { WeightType::DY, "DY" },
    { WeightType::TTbar, "TTbar" },
    { WeightType::Wjets, "Wjets" },
    { WeightType::BSM_to_SM, "BSM_to_SM" },
    { WeightType::TopPt, "TopPt" }
};

using WeightingMode = std::set<WeightType>;

inline WeightingMode operator|(const WeightingMode& a, const WeightingMode& b)
{
    WeightingMode c(a.begin(), a.end());
    c.insert(b.begin(), b.end());
    return c;
}

inline WeightingMode operator&(const WeightingMode& a, const WeightingMode& b)
{
    WeightingMode c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.end()));
    return c;
}

} // namespace analysis
} // namespace mc_corrections
