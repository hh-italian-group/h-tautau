/*! Jet PU Id weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class JetPuIdWeights : public IWeightProvider {
public:
    JetPuIdWeights(const std::string& file_eff, const std::string& file_sf, const std::string& file_mistag_eff,
                   const std::string& file_mistag_sf, const BTagger& _bTagger, Period _period);
    double GetEfficiency(std::shared_ptr<TH2F> hist, double pt, double eta,
                         UncertaintySource unc_source = UncertaintySource::None,
                         UncertaintyScale unc_scale = UncertaintyScale::Central) const;

    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;

    double GetWeight(EventInfo& eventInfo, UncertaintySource unc_source = UncertaintySource::None,
                     UncertaintyScale unc_scale = UncertaintyScale::Central) const;

private:
    std::shared_ptr<TH2F> eff_hist;
    std::shared_ptr<TH2F> sf_hist;
    std::shared_ptr<TH2F> eff_mistag_hist;
    std::shared_ptr<TH2F> sf_mistag_hist;
    BTagger bTagger;
    Period period;
};
} // namespace mc_corrections
} // namespace analysis
