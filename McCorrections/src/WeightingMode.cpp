/*! Weighting mode.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/McCorrections/include/WeightingMode.h"

namespace analysis {
namespace mc_corrections {

WeightingMode WeightingMode::operator|(const WeightingMode& b) const
{
    WeightingMode c;
    c.insert(this->begin(), this->end());
    c.insert(b.begin(), b.end());
    return c;
}

WeightingMode WeightingMode::operator&(const WeightingMode& b) const
{
    WeightingMode c;
    std::set_intersection(this->begin(), this->end(), b.begin(), b.end(), std::inserter(c, c.end()));
    return c;
}

WeightingMode WeightingMode::operator-(const WeightingMode& b) const
{
    WeightingMode c;
    std::set_difference(this->begin(), this->end(), b.begin(), b.end(), std::inserter(c, c.end()));
    return c;
}

} // namespace analysis
} // namespace mc_corrections
