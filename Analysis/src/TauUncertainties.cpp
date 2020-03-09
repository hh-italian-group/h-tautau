/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"
// #include "TGraphAsymmErrors.h"
// #include "AnalysisTools/Core/include/RootExt.h"
// #include <string>

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

std::shared_ptr<TH1F> TauESUncertainties::LoadWeight(const std::string& file_name)
{
    auto file = root_ext::OpenRootFile(file_name);
    return std::shared_ptr<TH1F>(root_ext::ReadCloneObject<TH1F>(*file, "tes", "", true));
}

TauESUncertainties::TauESUncertainties(std::string file_tes_low_pt, std::string file_tes_high_pt,
                                       std::string files_ele_faking_tau):
hist_tes_pt_low(LoadWeight(file_tes_low_pt)), hist_tes_pt_high(LoadWeight(file_tes_high_pt))
{

    // auto file_low = root_ext::OpenRootFile(file_tes_low_pt);
    // hist_tes_pt_low = root_ext::ReadCloneObject<TH1F>(*file_low, "tes");
    // auto file_high = root_ext::OpenRootFile(file_tes_high_pt);
    // hist_tes_pt_high = root_ext::ReadCloneObject<TH1F>(*file_high, "tes");

    // std::vector<int> dms_tes = {0, 1, 10, 11};
    // const float pt_low  = 34;  // average pT in Z -> tautau measurement (incl. in DM)
    // const float pt_high = 170; // average pT in W* -> taunu measurement (incl. in DM)
    // for (size_t dm = 0; dm < dms_tes.size(); ++dm){
    //     Int_t bin = hist_tes_pt_low->GetXaxis()->FindBin(dms_tes.at(dm));
    //     double tes  = hist_tes_pt_low->GetBinContent(bin);
    //
    //     double err = 0;
    //     if(pt > pt_high){
    //         Int_t bin_high = hist_tes_pt_high->GetXaxis()->FindBin(dms_tes.at(dm));
    //         err = hist_tes_pt_high->GetBinError(bin_high);
    //     }
    //     else if(pt > pt_low){
    //         Int_t bin_high = hist_tes_pt_high->GetXaxis()->FindBin(dms_tes.at(dm));
    //         double err_high = hist_tes_pt_high->GetBinError(bin_high);
    //
    //         double err_low  = hist_tes_pt_low->GetBinError(bin);
    //         err = err_low + (err_high - err_low)/(pt_high - pt_low)*(pt - pt_low);
    //     }
    //     else
    //         err = hist_tes_pt_low->GetBinError(bin);
    //
    //     tes_map[dms_tes.at(dm)] = PhysicalValue(tes, err);
    // }

    auto file_ele_faking_tau = root_ext::OpenRootFile(files_ele_faking_tau);
    auto hist_ele_faking_tau = root_ext::ReadObject<TGraphAsymmErrors>(*file_ele_faking_tau, "fes");

    std::vector<int> dms_fes = {0, 1};

    Int_t i = 0;

    for (int n_bin = 0; n_bin < 2; ++n_bin){
        for (int dm = 0; dm < dms_fes.size(); ++dm){
            bool is_bin_barrel = (i == 0 || i == 1 ) ? true : false;
            double y = hist_ele_faking_tau->GetY()[i];
            double y_error_up = hist_ele_faking_tau->GetErrorYhigh(i);
            double y_error_low = hist_ele_faking_tau->GetErrorYlow(i);
            fes[std::make_pair(dm, is_bin_barrel)] = y;
            fes_error[std::make_pair(dm, is_bin_barrel)] = std::make_pair(y_error_low, y_error_up);
            ++i;
        }
    }
}

double TauESUncertainties::GetCorrectionFactor(int decayMode, GenLeptonMatch genLeptonMatch,
                                               UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                               double eta) const
{
    if(genLeptonMatch == GenLeptonMatch::Tau) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorTrueTau(pt, decayMode, current_scale, genLeptonMatch);

    }
    else if(genLeptonMatch == GenLeptonMatch::Electron || genLeptonMatch == GenLeptonMatch::TauElectron) {
        UncertaintyScale current_scale = unc_source == UncertaintySource::EleFakingTauES
                                       ? scale : UncertaintyScale::Central;
                                       // return 1;
        return GetCorrectionFactorEleFakingTau(current_scale, eta, genLeptonMatch, decayMode);
    }
    else if(genLeptonMatch == GenLeptonMatch::Muon || genLeptonMatch == GenLeptonMatch::TauMuon){
        UncertaintyScale current_scale = unc_source == UncertaintySource::MuFakingTauES
                                       ? scale : UncertaintyScale::Central;
        return GetCorrectionFactorMuonFakingTau(current_scale);
    }
    else
        return 1.;

}

double TauESUncertainties::GetCorrectionFactorTrueTau(double pt, int decayMode, UncertaintyScale scale,
                                                      GenLeptonMatch genLeptonMatch) const
{
    std::vector<int> dms = {0, 1, 10, 11};
    const float pt_low  = 34;  // average pT in Z -> tautau measurement (incl. in DM)
    const float pt_high = 170; // average pT in W* -> taunu measurement (incl. in DM)
    if(std::find(dms.begin(), dms.end(), decayMode) != dms.end()){
        Int_t bin = hist_tes_pt_low->GetXaxis()->FindBin(decayMode);
        double tes  = hist_tes_pt_low->GetBinContent(bin);

        double err = 0;
        if(pt > pt_high){
            Int_t bin_high = hist_tes_pt_high->GetXaxis()->FindBin(decayMode);
            err = hist_tes_pt_high->GetBinError(bin_high);
        }
        else if(pt > pt_low){
            Int_t bin_high = hist_tes_pt_high->GetXaxis()->FindBin(decayMode);
            double err_high = hist_tes_pt_high->GetBinError(bin_high);

            double err_low  = hist_tes_pt_low->GetBinError(bin);
            err = err_low + (err_high - err_low)/(pt_high - pt_low)*(pt - pt_low);
        }
        else
            err = hist_tes_pt_low->GetBinError(bin);

        double tau_final_correction = tes + static_cast<int>(scale) * err; //* tes_error;
        return tau_final_correction;
    }
    return 1.;
}

double TauESUncertainties::GetCorrectionFactorEleFakingTau(UncertaintyScale scale, double eta,
                                                           GenLeptonMatch genLeptonMatch, int decayMode) const
{
    // if(std::find(dms.begin(), dms.end(), decayMode) != dms.end()){
    if(decayMode == 0 || decayMode == 1){
        bool is_barrel = std::abs(eta) < 1.5 ? true : false;

        auto errors = fes_error.at(std::make_pair(decayMode, is_barrel));
        auto error_fes  = static_cast<int>(scale) > 0 ? errors.second : errors.first;
        auto fes_value = fes.at(std::make_pair(decayMode, is_barrel));

        auto e_fake_rate_final_correction = fes_value + static_cast<int>(scale) * error_fes;
        return e_fake_rate_final_correction;
    }
    return 1.;
}


double TauESUncertainties::GetCorrectionFactorMuonFakingTau(UncertaintyScale scale) const
{
    return 1. + static_cast<int>(scale) * 0.01;
}
} // namespace analysis
