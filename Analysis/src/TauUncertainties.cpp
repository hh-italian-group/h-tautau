/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/TauUncertainties.h"
#include <TGraphAsymmErrors.h>
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {

std::map<int, StVariable> TauESUncertainties::LoadTauCorrections(const std::string& file_name)
{
    static const std::vector<int> available_decay_modes = { 0, 1, 10, 11 };
    auto file = root_ext::OpenRootFile(file_name);
    auto hist = std::shared_ptr<TH1F>(root_ext::ReadObject<TH1F>(*file, "tes"));
    std::map<int, StVariable> corr_map;
    for(int dm : available_decay_modes) {
        const int bin = hist->GetXaxis()->FindBin(dm);
        corr_map[dm] = StVariable(hist->GetBinContent(bin), hist->GetBinError(bin));
    }
    return corr_map;
}

std::map<std::pair<int, bool>, StVariable> TauESUncertainties::LoadElectronCorrections(const std::string& file_name)
{
    static const std::vector<int> available_decay_modes = { 0, 1 };
    auto file = root_ext::OpenRootFile(file_name);
    auto graph = std::shared_ptr<TGraphAsymmErrors>(root_ext::ReadObject<TGraphAsymmErrors>(*file, "fes"));
    std::map<std::pair<int, bool>, StVariable> corr_map;
    int n = 0;
    for(size_t region_id = 0; region_id < 2; ++region_id) {
        for(int dm : available_decay_modes) {
            const bool is_barrel = region_id == 0;
            corr_map[std::make_pair(dm, is_barrel)] = StVariable(graph->GetY()[n],
                                                                graph->GetErrorYhigh(n),
                                                                graph->GetErrorYlow(n));
            ++n;
        }
    }
    return corr_map;
}

bool TauESUncertainties::ApplyUncertaintyScale(int decayMode, GenLeptonMatch genLeptonMatch,
                                               UncertaintySource unc_source)
{
    static const std::map<GenLeptonMatch, std::map<UncertaintySource, int>> unc_source_map = {
        { GenLeptonMatch::Tau, { { UncertaintySource::TauES, -1 },
            { UncertaintySource::TauES_DM0, 0 }, { UncertaintySource::TauES_DM1, 1 },
            { UncertaintySource::TauES_DM10, 10 }, { UncertaintySource::TauES_DM11, 11 } }
        },
        { GenLeptonMatch::Electron, { { UncertaintySource::EleFakingTauES, -1 },
            { UncertaintySource::EleFakingTauES_DM0, 0 }, { UncertaintySource::EleFakingTauES_DM1, 0 } }
        },
        { GenLeptonMatch::TauElectron, { { UncertaintySource::EleFakingTauES, -1 },
            { UncertaintySource::EleFakingTauES_DM0, 0 }, { UncertaintySource::EleFakingTauES_DM1, 0 } }
        },
        { GenLeptonMatch::Muon, { { UncertaintySource::MuFakingTauES, -1 },  } },
        { GenLeptonMatch::TauMuon, { { UncertaintySource::MuFakingTauES, -1 },  } },
    };

    const auto gen_iter = unc_source_map.find(genLeptonMatch);
    if(gen_iter != unc_source_map.end()) {
        const auto source_iter = gen_iter->second.find(unc_source);
        if(source_iter != gen_iter->second.end()) {
            return source_iter->second < 0 || decayMode == source_iter->second;
        }
    }
    return false;
}

TauESUncertainties::TauESUncertainties(const std::string& file_tes_low_pt, const std::string& file_tes_high_pt,
                                       const std::string& file_ele_faking_tau) :
    tes_low_pt(LoadTauCorrections(file_tes_low_pt)), tes_high_pt(LoadTauCorrections(file_tes_high_pt)),
    fes(LoadElectronCorrections(file_ele_faking_tau))
{
}

double TauESUncertainties::GetCorrectionFactor(int decayMode, GenLeptonMatch genLeptonMatch,
                                               UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                               double eta, bool* same_as_central) const
{
    const UncertaintyScale current_scale = ApplyUncertaintyScale(decayMode, genLeptonMatch, unc_source)
                                         ? scale : UncertaintyScale::Central;
    if(same_as_central)
        *same_as_central = current_scale == UncertaintyScale::Central;
    if(genLeptonMatch == GenLeptonMatch::Tau) {
        return GetCorrectionFactorTrueTau(pt, decayMode, current_scale);

    }
    else if(genLeptonMatch == GenLeptonMatch::Electron || genLeptonMatch == GenLeptonMatch::TauElectron) {
        return GetCorrectionFactorEleFakingTau(eta, decayMode, current_scale);
    }
    else if(genLeptonMatch == GenLeptonMatch::Muon || genLeptonMatch == GenLeptonMatch::TauMuon){
        return GetCorrectionFactorMuonFakingTau(current_scale);
    }
    else {
        if(same_as_central)
            *same_as_central = true;
        return 1.;
    }
}

double TauESUncertainties::GetCorrectionFactorTrueTau(double pt, int decayMode, UncertaintyScale scale) const
{
    static constexpr double pt_low  = 34;  // average pT in Z -> tautau measurement (incl. in DM)
    static constexpr double pt_high = 170; // average pT in W* -> taunu measurement (incl. in DM)

    const auto low_pt_iter = tes_low_pt.find(decayMode);
    if(low_pt_iter == tes_low_pt.end())
        return 1.;
    double tes = low_pt_iter->second.value;
    const double err_low = low_pt_iter->second.error_up;
    double err = err_low;
    if(pt > pt_low) {
        const double err_high = tes_high_pt.at(decayMode).error_up;
        if(pt > pt_high)
            err = err_high;
        else
            err = err_low + (err_high - err_low)/(pt_high - pt_low)*(pt - pt_low);
    }

    return tes + static_cast<int>(scale) * err;
}

double TauESUncertainties::GetCorrectionFactorEleFakingTau(double eta, int decayMode, UncertaintyScale unc_scale) const
{
    const auto iter = fes.find(std::make_pair(decayMode, std::abs(eta) < 1.5));
    if(iter == fes.end())
        return 1.;
    const int scale = static_cast<int>(unc_scale);
    const double fes_value = iter->second.value;
    const double fes_error = scale > 0 ? iter->second.error_up : iter->second.error_low;
    return fes_value + scale * fes_error;
}

double TauESUncertainties::GetCorrectionFactorMuonFakingTau(UncertaintyScale scale) const
{
    return 1. + static_cast<int>(scale) * 0.01;
}

} // namespace analysis
