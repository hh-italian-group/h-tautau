/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class GenEventWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    
    virtual double Get(const Event& event) const override { return event.genEventWeight; }
    virtual double Get(const ntuple::ExpressEvent& event) const override { return event.genEventWeight; }
};
} // namespace mc_corrections
} // namespace analysis
