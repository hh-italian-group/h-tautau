/*! Jet PU Id weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

struct JetInfo {
    double pt, eta;
    double eff, SF;
    bool jetPuIdOutcome;

    JetInfo(const JetCandidate& jet);
};

class JetPuIdWeights : public IWeightProvider {
public:
    JetPuIdWeights(const std::string& file_eff, const std::string& file_sf, const BTagger& _bTagger, Period _period);
    double GetEfficiency(std::shared_ptr<TH2F> hist, double pt, double eta) const;
    void InitializeEff(JetInfo& jetInfo) const;

    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;

private:
    static double GetJetPuIdWeight(const std::vector<JetInfo>& jetInfos);


private:
    std::shared_ptr<TH2F> eff_hist;
    std::shared_ptr<TH2F> sf_hist;
    BTagger bTagger;
    Period period;
};
} // namespace mc_corrections
} // namespace analysis
