/*! Various event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include <boost/filesystem.hpp>
#include "PileUpWeight.h"
#include "LeptonWeights.h"
#include "BTagWeight.h"
#include "TopPtWeight.h"
#include "TauIdWeight.h"
#include "WeightingMode.h"
#include "h-tautau/Analysis/include/EventInfo.h"


namespace analysis {
namespace mc_corrections {

class EventWeights {
public:
    using ProviderPtr = std::shared_ptr<IWeightProvider>;
    using ProviderMap = std::map<WeightType, ProviderPtr>;

    EventWeights(Period period, JetOrdering jet_ordering, DiscriminatorWP btag_wp, const WeightingMode& mode = {})
    {
        if(period == Period::Run2016) {
            if(mode.empty() || mode.count(WeightType::PileUp))
                providers[WeightType::PileUp] = std::make_shared<PileUpWeight>(
                            FullName("2016/pileup_weight_600bins_Moriond17.root"), "pileup_weight", 60, 0);
            if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
                providers[WeightType::LeptonTrigIdIso] = std::make_shared<LeptonWeights>(
                            FullLeptonName("Electron/Run2016BtoH/Electron_IdIso_IsoLt0p15_eff.root"),
                            FullLeptonName("Electron/Run2016BtoH/Electron_Ele25WPTight_eff.root"),
                            "",
                            FullLeptonName("Muon/Run2016BtoH/Muon_IdIso_IsoLt0p2_2016BtoH_eff_update1407.root"),
                            FullLeptonName("Muon/Run2016BtoH/Muon_Mu22OR_eta2p1_eff.root"), "",
                            FullName("2016/Tau/fitresults_tt_moriond2017.json"), "",period, DiscriminatorWP::Medium);
            if(mode.empty() || mode.count(WeightType::BTag))
                providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullName("2016/btag/bTagEfficiencies_Moriond17.root"), FullName("2016/btag/CSVv2_Moriond17_B_H.csv"),
                            period, jet_ordering, btag_wp);
            if(mode.empty() || mode.count(WeightType::TopPt))
                providers[WeightType::TopPt] = std::make_shared<TopPtWeight>(0.0615, 0.0005);
        }

        else if(period == Period::Run2017) {
            if(mode.empty() || mode.count(WeightType::PileUp))
                providers[WeightType::PileUp] = std::make_shared<PileUpWeight>(
                            FullName("2017/pileup_weight_200bins.root"), "pileup_weight", 130, 0);
            if(mode.empty() || mode.count(WeightType::BTag)){
                if(jet_ordering == JetOrdering::DeepCSV)
                    providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullName("2017/btag/BTagEfficiency_deep_csv_pu_id_full.root"), FullName("2017/btag/DeepCSV_94XSF_V3_B_F.csv"),
                            period, jet_ordering, btag_wp);
                else if(jet_ordering == JetOrdering::CSV)
                    providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullName("2017/btag/BTagEfficiency_csv_pu_id_full.root"), FullName("2017/btag/CSVv2_94XSF_V2_B_F.csv"),
                            period, jet_ordering, btag_wp);
                else if(jet_ordering == JetOrdering::DeepFlavour)
                    providers[WeightType::BTag] = std::make_shared<BTagWeight>(
                            FullName("2017/btag/BTagEfficiency_deep_Flavour_pu_id_full.root"), FullName("2017/btag/DeepFlavour_94XSF_V1_B_F.csv"),
                            period, jet_ordering, btag_wp);
                else
                   throw exception("Jet_Ordering %1% is not supported.") % jet_ordering;
            }

                        if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
                        providers[WeightType::LeptonTrigIdIso] = std::make_shared<LeptonWeights>(
                                    FullLeptonName("Electron/Run2017/Electron_IdIso_IsoLt0.10_eff_RerecoFall17.root"),
                                    FullLeptonName("Electron/Run2017/Electron_Ele32orEle35.root"),
                                    FullLeptonName("Electron/Run2017/Electron_EleTau_Ele24.root"),
                                    FullLeptonName("Muon/Run2017/Muon_IdIso_IsoLt0.15_eff_RerecoFall17.root"),
                                    FullLeptonName("Muon/Run2017/Muon_IsoMu24orIsoMu27.root"),
                                    FullLeptonName("Muon/Run2017/Muon_MuTau_IsoMu20.root"),
                                    FullName("2017/Tau/tauTriggerEfficiencies2017_New.root"),
                                    FullName("2017/Tau/tauTariggerEfficiencies.root"),
                                    period, DiscriminatorWP::Medium);
                                    // POG SFs
                                    // FullName("2017/Electron/EleIdSFPOG.root"),
                                    // FullName("2017/Electron/EleIsoSFPOG.root"),
                                    // FullName("2017/Muon/MuonIdSFPOG.root"),
                                    // FullName("2017/Muon/MuonIsoSFPOG.root"),
                                    // FullName("2017/Muon/MuonTriggerSFPOG.root"),
                                    // FullName("Tau/fitresults_tt_moriond2017.json"),
                    if(mode.empty() || mode.count(WeightType::TopPt))
                        providers[WeightType::TopPt] = std::make_shared<TopPtWeight>(0.0615, 0.0005);
        }
        else {
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
    static std::string FullName(const std::string& fileName, const std::string& path)
    {
        const std::string name = path + "/" + fileName;
        if (!boost::filesystem::exists(name))
            throw exception ("file '%1%' not found") % name;
        return name;
    }

    static std::string FullName(const std::string& fileName)
    {
        static const std::string path = "h-tautau/McCorrections/data";
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
