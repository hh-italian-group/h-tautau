/*! Various event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "PileUpWeight.h"
#include "LeptonWeights.h"
#include "BTagWeight.h"
#include "TopPtWeight.h"
#include "TauIdWeight.h"
#include "WeightingMode.h"

namespace analysis {
namespace mc_corrections {

class EventWeights {
public:
    using ProviderPtr = std::shared_ptr<IWeightProvider>;
    using ProviderMap = std::map<WeightType, ProviderPtr>;

    EventWeights(Period period, DiscriminatorWP btag_wp, const WeightingMode& mode = {})
    {
        if(period == Period::Run2015) {
            if(mode.empty() || mode.count(WeightType::PileUp))
                providers[WeightType::PileUp] = std::make_shared<PileUpWeight>(
                            FullName("reWeight_Fall.root"), "lumiWeights", 60, 0);
            if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
                providers[WeightType::LeptonTrigIdIso] = std::make_shared<LeptonWeights>(
                            FullLeptonName("Electron/Electron_IdIso0p10_eff.root"),
                            FullLeptonName("Electron/Electron_SingleEle_eff.root"),
                            FullLeptonName("Muon/Muon_IdIso0p1_fall15.root"),
                            FullLeptonName("Muon/Muon_IsoMu18_fall15.root"));
            if(mode.empty() || mode.count(WeightType::BTag))
                providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullName("bTagEff_Loose.root"), FullName("CSVv2.csv"), btag_wp);
        }
        else if(period == Period::Run2016) {
            if(mode.empty() || mode.count(WeightType::PileUp))
                providers[WeightType::PileUp] = std::make_shared<PileUpWeight>(
                            FullName("pileup_weight_600bins_Moriond17.root"), "pileup_weight", 60, 0);
            if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
                providers[WeightType::LeptonTrigIdIso] = std::make_shared<LeptonWeights>(
                            FullLeptonName("Electron/Run2016BtoH/Electron_IdIso_IsoLt0p15_eff.root"),
                            FullLeptonName("Electron/Run2016BtoH/Electron_Ele25WPTight_eff.root"),
                            FullLeptonName("Muon/Run2016BtoH/Muon_IdIso_IsoLt0p2_2016BtoH_eff_update1407.root"),
                            FullLeptonName("Muon/Run2016BtoH/Muon_Mu22OR_eta2p1_eff.root"));
            if(mode.empty() || mode.count(WeightType::BTag))
                providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullBtagName("bTagEfficiencies_Moriond17.root"), FullBtagName("CSVv2_Moriond17_B_H.csv"),
                            btag_wp);
            if(mode.empty() || mode.count(WeightType::TopPt))
                providers[WeightType::TopPt] = std::make_shared<TopPtWeight>(0.0615, 0.0005);
            if(mode.empty() || mode.count(WeightType::TauId))
                providers[WeightType::TauId] = std::make_shared<TauIdWeight>(
                            FullName("Tau/fitresults_tt_moriond2017.json"), DiscriminatorWP::Medium);
        } else {
            throw exception("Period %1% is not supported.") % period;
        }
    }

    ProviderPtr GetProvider(WeightType weightType) const
    {
        if(!providers.count(weightType))
            throw exception("Weight provider not found for %1% weight.") % weightType;
        return providers.at(weightType);
    }

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
        return GetProvider(weightType)->Get(event);
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
    static std::string FullName(const std::string& fileName, const std::string& path)
    {
        return path + "/" + fileName;
    }

    static std::string FullName(const std::string& fileName)
    {
        static const std::string path = "h-tautau/McCorrections/data";
        return FullName(fileName, path);
    }

    static std::string FullBtagName(const std::string& fileName)
    {
        static const std::string path = "h-tautau/McCorrections/data/btag";
        return FullName(fileName, path);
    }

    static std::string FullLeptonName(const std::string& fileName)
    {
        static const std::string path = "HTT-utilities/LepEffInterface/data";
        return FullName(fileName, path);
    }

protected:
    ProviderMap providers;
};

} // namespace mc_corrections
} // namespace analysis
