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

    DiscriminatorIdResults() : results(0) {}
    DiscriminatorIdResults(BitsContainer _results) : results(_results) {}

    bool Passed(DiscriminatorWP wp) const
    {
        const unsigned bit_index = static_cast<unsigned>(wp);
        if(bit_index > MaxNumberOfWorkingPoints)
            throw exception("Discriminator WP = '{}' is not supported.") % wp;

        const BitsContainer mask = static_cast<BitsContainer>(BitsContainer(1) << bit_index);
        return (results & mask) != BitsContainer(0);
    }
    bool Failed(DiscriminatorWP wp) const { return !Passed(wp); }

    void SetResult(DiscriminatorWP wp, bool result)
    {
        const unsigned bit_index = static_cast<unsigned>(wp);
        if(bit_index > MaxNumberOfWorkingPoints)
            throw exception("Discriminator WP = '{}' is not supported.") % wp;
        const BitsContainer mask = static_cast<BitsContainer>(BitsContainer(1) << bit_index);
        results = (results & ~mask) | static_cast<BitsContainer>(BitsContainer(result) << bit_index);
    }

    BitsContainer GetResultBits() const { return results; }

private:
    BitsContainer results;
};

} // namespace analysis
