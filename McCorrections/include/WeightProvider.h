/*! Weight provider interface.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"

namespace analysis {
namespace mc_corrections {

class IWeightProvider {
public:
    virtual ~IWeightProvider() {}
    virtual double Get(const ntuple::Event& event) const = 0;
    virtual double Get(const ntuple::ExpressEvent& event) const = 0;
};

} // namespace mc_corrections
} // namespace analysis
