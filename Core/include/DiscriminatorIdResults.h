/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include "AnalysisTypes.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

struct DiscriminatorIdResults {
    using BitsContainer = uint16_t;
    static constexpr size_t MaxNumberOfWorkingPoints = std::numeric_limits<BitsContainer>::digits;

    DiscriminatorIdResults();
    DiscriminatorIdResults(BitsContainer _results);

    bool Passed(DiscriminatorWP wp) const;
    bool Failed(DiscriminatorWP wp) const;

    void SetResult(DiscriminatorWP wp, bool result);

    BitsContainer GetResultBits() const;

private:
    BitsContainer results;
};

} // namespace analysis
