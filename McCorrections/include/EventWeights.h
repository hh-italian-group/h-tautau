/*! Various event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/McCorrections/include/WeightingMode.h"
#include "h-tautau/McCorrections/include/WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class EventWeights {
public:
    using ProviderPtr = std::shared_ptr<IWeightProvider>;
    using ProviderMap = std::map<WeightType, ProviderPtr>;

    EventWeights(Period period, JetOrdering jet_ordering, DiscriminatorWP btag_wp, const WeightingMode& mode = {});
    ProviderPtr GetProvider(WeightType weightType) const;

    template<typename Provider>
    std::shared_ptr<Provider> GetProviderT(WeightType weightType) const
    {
        auto provider = GetProvider(weightType);
        auto casted_provider = std::dynamic_pointer_cast<Provider>(provider);
        if(!casted_provider)
            throw exception("Can't cast provider for weight %1% to type %2%.") % weightType % typeid(Provider).name();
        return casted_provider;
    }

    template<typename Event>
    double GetWeight(const Event& event, WeightType weightType) const
    {
        double weight = GetProvider(weightType)->Get(event);
        if (std::isnan(weight) || std::abs(weight) == std::numeric_limits<double>::infinity())
            throw exception("%1% weights is nan or infinity for event %2%.") % weightType % EventIdentifier(event);
		return weight;
    }

    template<typename Event>
    double GetTotalWeight(const Event& event, const WeightingMode& weightingMode) const
    {
        double weight = 1.;
        for(WeightType weightType : weightingMode)
            weight *= GetWeight(event, weightType);
        return weight;
    }

protected:
    static std::string FullName(const std::string& fileName, const std::string& path);
    static std::string FullName(const std::string& fileName);
    static std::string FullLeptonName(const std::string& fileName);

protected:
    ProviderMap providers;
};

} // namespace mc_corrections
} // namespace analysis
