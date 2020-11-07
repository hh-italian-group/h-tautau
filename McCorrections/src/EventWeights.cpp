/*! Various event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/McCorrections/include/EventWeights.h"

#include <boost/filesystem.hpp>

#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/McCorrections/include/PileUpWeight.h"
#include "h-tautau/McCorrections/include/LeptonWeights.h"
#include "h-tautau/McCorrections/include/BTagWeight.h"
#include "h-tautau/McCorrections/include/TopPtWeight.h"
#include "h-tautau/McCorrections/include/GenEventWeight.h"
#include "h-tautau/McCorrections/include/JetPuIdWeights.h"

namespace analysis {
namespace mc_corrections {

EventWeights::EventWeights(Period period, const BTagger& bTagger, const WeightingMode& mode)
{
    static constexpr DiscriminatorWP default_btag_wp = DiscriminatorWP::Medium;

    const BTaggerKind base_tagger = bTagger.GetBaseTagger();

    if(period == Period::Run2016) {
        if(mode.empty() || mode.count(WeightType::PileUp))
            CreateProvider<PileUpWeightEx>(WeightType::PileUp,
                        FullName("2016/Pileup_Data2016.root"), FullName("2016/Pileup_Data2016_Up.root"),
                        FullName("2016/Pileup_Data2016_Down.root"),
                        FullName("2016/pu_mc_distr_per_sample_100_100_2016.root"),
                        FullName("2016/pileup_groups_2016.txt"), 100, 0);
        if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
            CreateProvider<LeptonWeights>(WeightType::LeptonTrigIdIso,
                        FullLeptonName("Electron/Run2016_legacy/Electron_Run2016_legacy_IdIso.root"),
                        FullLeptonName("Electron/Run2016_legacy/Electron_Run2016_legacy_Ele25.root"),
                        "",
                        FullLeptonName("Muon/Run2016_legacy/Muon_Run2016_legacy_IdIso.root"),
                        FullLeptonName("Muon/Run2016_legacy/Muon_Run2016_legacy_IsoMu22.root"),
                        FullLeptonName("Muon/Run2016BtoH/Muon_Mu19leg_2016BtoH_eff.root"),
                        FullTriggerName("2016_tauTriggerEff_DeepTau2017v2p1.root"),
                        period, false);
        if(mode.empty() || mode.count(WeightType::BTag)){
            if(base_tagger == BTaggerKind::CSV)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2016/btag/bTagEfficiencies_Moriond17.root"),
                        FullName("2016/btag/CSVv2_Moriond17_B_H.csv"),
                        bTagger, default_btag_wp);
            else if(base_tagger == BTaggerKind::DeepCSV)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2016/btag/b_eff_HH_DeepCSV_2016.root"),
                        FullName("2016/btag/DeepCSV_2016LegacySF_WP_V1.csv"),
                        bTagger, default_btag_wp);
            else if(base_tagger == BTaggerKind::DeepFlavour)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2016/btag/b_eff_HH_DeepFlavour_2016.root"),
                        FullName("2016/btag/DeepJet_2016LegacySF_WP_V1.csv"),
                        bTagger, default_btag_wp);
            else
                throw exception("EventWeights: b tagger %1% is not supported.") % base_tagger;
        }
        if(mode.empty() || mode.count(WeightType::JetPuIdWeights))
            CreateProvider<JetPuIdWeights>(WeightType::JetPuIdWeights,
                        FullName("jet_pu_id/effcyPUID_81Xtraining.root"),
                        FullName("jet_pu_id/scalefactorsPUID_81Xtraining.root"),
                        bTagger, period);
        if(mode.empty() || mode.count(WeightType::TopPt))
            CreateProvider<TopPtWeight>(WeightType::TopPt, 0.0615, 0.0005);
    }

    else if(period == Period::Run2017) {
        if(mode.empty() || mode.count(WeightType::PileUp))
            CreateProvider<PileUpWeightEx>(WeightType::PileUp,
                        FullName("2017/Pileup_Data2017.root"), FullName("2017/Pileup_Data2017_Up.root"),
                        FullName("2017/Pileup_Data2017_Down.root"),
                        FullName("2017/pu_mc_distr_per_sample_100_100_2017.root"),
                        FullName("2017/pileup_groups.txt"), 100, 0);
        if(mode.empty() || mode.count(WeightType::BTag)){
            if(base_tagger == BTaggerKind::DeepCSV)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2017/btag/b_eff_HH_DeepCSV_2017.root"),
                        FullName("2017/btag/DeepCSV_94XSF_WP_V4_B_F.csv"),
                        bTagger, default_btag_wp);
            else if(base_tagger == BTaggerKind::CSV)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2017/btag/BTagEfficiency_csv_pu_id_full.root"),
                        FullName("2017/btag/CSVv2_94XSF_V2_B_F.csv"),
                        bTagger, default_btag_wp);
            else if(base_tagger == BTaggerKind::DeepFlavour)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2017/btag/b_eff_HH_DeepFlavour_2017.root"),
                        FullName("2017/btag/DeepFlavour_94XSF_WP_V3_B_F.csv"),
                        bTagger, default_btag_wp);
            else
               throw exception("EventWeights: b tagger %1% is not supported.") % base_tagger;
        }
        if(mode.empty() || mode.count(WeightType::JetPuIdWeights))
            CreateProvider<JetPuIdWeights>(WeightType::JetPuIdWeights,
                        FullName("jet_pu_id/effcyPUID_81Xtraining.root"),
                        FullName("jet_pu_id/scalefactorsPUID_81Xtraining.root"),
                        bTagger, period);
        if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
            CreateProvider<LeptonWeights>(WeightType::LeptonTrigIdIso,
                        FullLeptonName("Electron/Run2017/Electron_Run2017_IdIso.root"),
                        FullLeptonName("Electron/Run2017/Electron_Ele32orEle35.root"),
                        FullLeptonName("Electron/Run2017/Electron_EleTau_Ele24.root"),
                        FullLeptonName("Muon/Run2017/Muon_IdIso_IsoLt0.15_eff_RerecoFall17.root"),
                        FullLeptonName("Muon/Run2017/Muon_IsoMu24orIsoMu27.root"),
                        FullLeptonName("Muon/Run2017/Muon_MuTau_IsoMu20.root"),
                        FullTriggerName("2017_tauTriggerEff_DeepTau2017v2p1.root"),
                        period, false);
        if(mode.empty() || mode.count(WeightType::TopPt))
            CreateProvider<TopPtWeight>(WeightType::TopPt, 0.0615, 0.0005);
    }
    else if(period == Period::Run2018) {
        if(mode.empty() || mode.count(WeightType::PileUp))
            CreateProvider<PileUpWeightEx>(WeightType::PileUp,
                        FullName("2018/Pileup_Data2018.root"), FullName("2018/Pileup_Data2018_Up.root"),
                        FullName("2018/Pileup_Data2018_Down.root"),
                        FullName("2018/pu_mc_distr_per_sample_100_100_2018.root"),
                        FullName("2018/pileup_groups_2018.txt"), 100, 0);
        if(mode.empty() || mode.count(WeightType::BTag)){
            if(base_tagger == BTaggerKind::DeepCSV)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2018/btag/b_eff_HH_DeepCSV_2018.root"),
                        FullName("2018/btag/DeepCSV_102XSF_WP_V1.csv"),
                        bTagger, default_btag_wp);
            else if(base_tagger == BTaggerKind::DeepFlavour)
                CreateProvider<BTagWeight>(WeightType::BTag,
                        FullName("2018/btag/b_eff_HH_DeepFlavour_2018.root"),
                        FullName("2018/btag/DeepJet_102XSF_V2.csv"),
                        bTagger, default_btag_wp);
            else
               throw exception("EventWeights: b tagger %1% is not supported.") % base_tagger;
        }
        if(mode.empty() || mode.count(WeightType::JetPuIdWeights))
            CreateProvider<JetPuIdWeights>(WeightType::JetPuIdWeights,
                        FullName("jet_pu_id/effcyPUID_81Xtraining.root"),
                        FullName("jet_pu_id/scalefactorsPUID_81Xtraining.root"),
                        bTagger, period);
        if(mode.empty() || mode.count(WeightType::TopPt))
            CreateProvider<TopPtWeight>(WeightType::TopPt, 0.0615, 0.0005);
        if(mode.empty() || mode.count(WeightType::LeptonTrigIdIso))
            CreateProvider<LeptonWeights>(WeightType::LeptonTrigIdIso,
                        FullLeptonName("Electron/Run2018/Electron_Run2018_IdIso.root"),
                        FullLeptonName("Electron/Run2018/Electron_Run2018_Ele32orEle35.root"),
                        FullLeptonName("Electron/Run2018/Electron_Run2018_Ele24.root"),
                        FullLeptonName("Muon/Run2018/Muon_Run2018_IdIso.root"),
                        FullLeptonName("Muon/Run2018/Muon_Run2018_IsoMu24orIsoMu27.root"),
                        FullLeptonName("Muon/Run2018/Muon_Run2018_IsoMu20.root"),
                        FullTriggerName("2018_tauTriggerEff_DeepTau2017v2p1.root"),
                        period, false);
    }
    else {
        throw exception("Period %1% is not supported (EventWeights).") % period;
    }
    if(mode.empty() || mode.count(WeightType::GenEventWeight))
        CreateProvider<GenEventWeight>(WeightType::GenEventWeight);
}

EventWeights::ProviderPtr EventWeights::GetProvider(WeightType weightType) const
{
    if(!providers.count(weightType))
        throw exception("Weight provider not found for %1% weight.") % weightType;
    return providers.at(weightType);
}

double EventWeights::GetWeight(EventInfo& event, WeightType weightType) const
{
    double weight = GetProvider(weightType)->Get(event);
    if (std::isnan(weight) || std::abs(weight) == std::numeric_limits<double>::infinity())
        throw exception("%1% weights is nan or infinity for event %2%.") % weightType % EventIdentifier(*event);
    return weight;
}


double EventWeights::GetTotalWeight(EventInfo& event, const WeightingMode& weightingMode) const
{
    double weight = 1.;
    for(WeightType weightType : weightingMode)
        weight *= GetWeight(event, weightType);
    return weight;
}

double EventWeights::GetWeight(const ntuple::ExpressEvent& event, WeightType weightType) const
{
    double weight = GetProvider(weightType)->Get(event);
    if (std::isnan(weight) || std::abs(weight) == std::numeric_limits<double>::infinity())
        throw exception("%1% weights is nan or infinity for event %2%.") % weightType % EventIdentifier(event);
    return weight;
}


double EventWeights::GetTotalWeight(const ntuple::ExpressEvent& event, const WeightingMode& weightingMode) const
{
    double weight = 1.;
    for(WeightType weightType : weightingMode)
        weight *= GetWeight(event, weightType);
    return weight;
}

std::string EventWeights::FullName(const std::string& fileName, const std::string& path)
{
    const std::string name = path + "/" + fileName;
    if (!boost::filesystem::exists(name))
        throw exception ("file '%1%' not found") % name;
    return name;
}

std::string EventWeights::FullName(const std::string& fileName)
{
    static const std::string path = "h-tautau/McCorrections/data";
    return FullName(fileName, path);
}

std::string EventWeights::FullLeptonName(const std::string& fileName)
{
    static const std::string path = "HTT-utilities/LepEffInterface/data";
    return FullName(fileName, path);
}

std::string EventWeights::FullTriggerName(const std::string& fileName)
{
    static const std::string path = "TauAnalysisTools/TauTriggerSFs/data/";
    return FullName(fileName, path);
}


} // namespace mc_corrections
} // namespace analysis
