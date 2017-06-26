/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

namespace detail {
class LeptonScaleFactors {
public:
    using SF = htt_utilities::ScaleFactor;

    LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInput) :
        idIso(new SF()), trigger(new SF())
    {
        idIso->init_ScaleFactor(idIsoInput);
        trigger->init_ScaleFactor(triggerInput);
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const
    {
        return idIso->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4) const
    {
        return trigger->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTotalSF(const LorentzVector& p4) const { return GetIsoSF(p4) * GetTriggerSF(p4); }

private:
    std::shared_ptr<SF> idIso, trigger;
};

} // namespace detail

class LeptonWeights : public IWeightProvider {
public:
    using Event = ntuple::Event;

    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_triggerInput,
                  const std::string& muon_idIsoInput, const std::string& muon_triggerInput) :
        electronSF(electron_idIsoInput, electron_triggerInput),
        muonSF(muon_idIsoInput, muon_triggerInput)
    {
    }

    double GetIdIsoWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::ETau) return electronSF.GetIdIsoSF(event.p4_1);
        if(channel == Channel::MuTau) return muonSF.GetIdIsoSF(event.p4_1);
        return 1.0;
    }

    double GetTriggerWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::ETau) return electronSF.GetTriggerSF(event.p4_1);
        if(channel == Channel::MuTau) return muonSF.GetTriggerSF(event.p4_1);
        return 1.0;
    }

    virtual double Get(const Event& event) const override { return GetIdIsoWeight(event) * GetTriggerWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    detail::LeptonScaleFactors electronSF, muonSF;
};

} // namespace mc_corrections
} // namespace analysis
