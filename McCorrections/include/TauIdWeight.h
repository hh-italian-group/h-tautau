/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int decay_mode,
                              DiscriminatorWP anti_ele_wp, DiscriminatorWP anti_mu_wp,
                              DiscriminatorWP iso_wp) const = 0;
    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP iso_wp, GenMatch /*gen_match*/) const = 0;
    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match,
                                            DiscriminatorWP anti_mu_wp) const = 0;
    virtual double GetEleMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match,
                                           DiscriminatorWP anti_ele_wp) const = 0;
    virtual ~TauIdWeight() {}
};

class TauIdWeight2016 : public TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/, int /*decay_mode*/,
                              DiscriminatorWP /*anti_ele_wp*/, DiscriminatorWP /*anti_mu_wp*/,
                              DiscriminatorWP /*iso_wp*/) const override;
    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP /*iso_wp*/, GenMatch /*gen_match*/) const override;
    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/,
                                            DiscriminatorWP /*anti_mu_wp*/) const override;
    virtual double GetEleMissIdUncertainty(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/,
                                           DiscriminatorWP /*anti_ele_wp*/) const override;
};

class TauIdWeight2017 : public TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int /*decay_mode*/,
                              DiscriminatorWP anti_ele_wp, DiscriminatorWP anti_mu_wp,
                              DiscriminatorWP iso_wp) const override;
    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP iso_wp, GenMatch /*gen_match*/) const override;
    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match,
                                            DiscriminatorWP anti_mu_wp) const override;
    virtual double GetEleMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match,
                                           DiscriminatorWP anti_ele_wp) const override;

public:
    PhysicalValue getMuonMissId(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP iso_wp) const;
    PhysicalValue getEleMissId(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP iso_wp) const;
    PhysicalValue getTauIso(DiscriminatorWP iso_wp, GenMatch /*gen_match*/) const;
    PhysicalValue tauIdForDM(GenMatch gen_match, int decay_mode) const;
};
} // namespace mc_corrections
} // namespace analysis
