/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/TauIdResults.h"
#include "h-tautau/Core/include/DiscriminatorIdResults.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

    DiscriminatorIdResults::DiscriminatorIdResults() : results(0) {}
    DiscriminatorIdResults::DiscriminatorIdResults(BitsContainer _results) : results(_results) {}

    bool DiscriminatorIdResults::Passed(DiscriminatorWP wp) const
    {
        const unsigned bit_index = static_cast<unsigned>(wp);
        if(bit_index > MaxNumberOfWorkingPoints)
            throw exception("Discriminator WP = '{}' is not supported.") % wp;

        const BitsContainer mask = static_cast<BitsContainer>(BitsContainer(1) << bit_index);
        return (results & mask) != BitsContainer(0);
    }
    bool DiscriminatorIdResults::Failed(DiscriminatorWP wp) const { return !Passed(wp); }

    void DiscriminatorIdResults::SetResult(DiscriminatorWP wp, bool result)
    {
        const unsigned bit_index = static_cast<unsigned>(wp);
        if(bit_index > MaxNumberOfWorkingPoints)
            throw exception("Discriminator WP = '{}' is not supported.") % wp;
        const BitsContainer mask = static_cast<BitsContainer>(BitsContainer(1) << bit_index);
        results = (results & ~mask) | static_cast<BitsContainer>(BitsContainer(result) << bit_index);
    }

    DiscriminatorIdResults::BitsContainer DiscriminatorIdResults::GetResultBits() const { return results; }

} // namespace analysis
