/*! Definition of weighting modes for H and HH analyses.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <set>
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {
namespace mc_corrections {

using ::analysis::operator<<;
using ::analysis::operator>>;

enum class WeightType {
    PileUp = 0, LeptonTrigIdIso = 1, BTag = 2, DY = 3, TTbar = 4, Wjets = 5, BSM_to_SM = 6, TopPt = 7, TauId = 8,
    GenEventWeight = 9, JetPuIdWeights = 10
};
ENUM_NAMES(WeightType) = {
    { WeightType::PileUp, "PileUp" },
    { WeightType::LeptonTrigIdIso, "LeptonTrigIdIso" },
    { WeightType::BTag, "BTag" },
    { WeightType::DY, "DY" },
    { WeightType::TTbar, "TTbar" },
    { WeightType::Wjets, "Wjets" },
    { WeightType::BSM_to_SM, "BSM_to_SM" },
    { WeightType::TopPt, "TopPt" },
    { WeightType::TauId, "TauId" },
    { WeightType::GenEventWeight, "GenEventWeight" },
    { WeightType::JetPuIdWeights, "JetPuIdWeights" }
};

class WeightingMode : public std::set<WeightType> {
public:
    using Base = std::set<WeightType>;

    WeightingMode() = default;
    WeightingMode(const WeightingMode&) = default;
    WeightingMode(WeightingMode&&) = default;
    WeightingMode& operator=(const WeightingMode&) = default;

    template<typename ...Args>
    WeightingMode(Args&&... args) { Initialize(std::forward<Args>(args)...); }

    WeightingMode operator|(const WeightingMode& b) const;
    WeightingMode operator&(const WeightingMode& b) const;
    WeightingMode operator-(const WeightingMode& b) const;

private:
    void Initialize() {}

    template<typename ...Args>
    void Initialize(WeightType value, Args&&... other_values)
    {
        insert(value);
        Initialize(std::forward<Args>(other_values)...);
    }
};

} // namespace analysis
} // namespace mc_corrections
