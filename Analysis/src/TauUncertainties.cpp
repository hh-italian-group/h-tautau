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

TauESUncertainties::TauESUncertainties(std::string file_tes_low_pt, std::string file_tes_high_pt,
                                       std::string files_ele_faking_tau, TauIdDiscriminator _tau_vs_e_discr):
    file_low(root_ext::OpenRootFile(file_tes_low_pt)),
    hist_tes_pt_low(root_ext::ReadCloneObject<TH1F>(*file_low, "tes")), //ReadCloneObject
    file_high(root_ext::OpenRootFile(file_tes_high_pt)), hist_tes_pt_high(root_ext::ReadObject<TH1F>(*file_high, "tes")),
    file_ele_faking_tau(root_ext::OpenRootFile(files_ele_faking_tau)),
    hist_ele_faking_tau(root_ext::ReadObject<TGraphAsymmErrors>(*file_ele_faking_tau, "fes")),
    tau_vs_e_discr(_tau_vs_e_discr)

{
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
        return GetCorrectionFactorEleFakingTau(current_scale, eta, genLeptonMatch, tau_vs_e_discr, decayMode);
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
    if(genLeptonMatch == GenLeptonMatch::Tau && (std::find(dms.begin(), dms.end(), decayMode) != dms.end())){
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

// double TauESUncertainties::GetCorrectionFactorEleFakingTau(UncertaintyScale scale, double eta,
//                                                            GenLeptonMatch genLeptonMatch,
//                                                            TauIdDiscriminator tauVSeDiscriminator, int decayMode) const
// {
//     std::vector<int> dms = {0, 1};
//     std::vector<std::string> regions = {"barrel", "endcap"};
//     std::map<std::pair<int, std::string>, double> fes; //mettere invece di string fare booleano
//     std::map<std::pair<int, std::string>, std::pair<double, double>> fes_error;
//
//     Int_t i = 0;
//     bool is_barrel = (i == 0 || i == 1 ) ? true : false;
//     for(const auto& region: regions){
//         for(const auto& dm: dms){
//             double y = hist_ele_faking_tau->GetY()[i];
//             double y_error_up = hist_ele_faking_tau->GetErrorYhigh(i);
//             double y_error_low = hist_ele_faking_tau->GetErrorYlow(i);
//             fes[std::make_pair(dm, region)] = y;
//             fes_error[std::make_pair(dm, region)] = std::make_pair(y_error_low, y_error_up);
//             ++i;
//         }
//     }
//
//     if((genLeptonMatch == GenLeptonMatch::Electron ||genLeptonMatch == GenLeptonMatch::TauElectron) &&
//         tauVSeDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSe &&
//         (std::find(dms.begin(), dms.end(), decayMode) != dms.end())){
//
//         std::string reg = "";
//         if(std::abs(eta) < 1.5)
//             reg = "barrel";
//         else
//             reg = "endcap";
//
//         auto errors = fes_error[std::make_pair(decayMode, reg)];
//         auto error_fes  = static_cast<int>(scale) > 0 ? errors.second : errors.first;
//         auto fes_value = fes[std::make_pair(decayMode, reg)];
//
//         auto e_fake_rate_final_correction = fes_value + static_cast<int>(scale) * error_fes;
//         return e_fake_rate_final_correction;
//     }
//     return 1.;
// }

double TauESUncertainties::GetCorrectionFactorEleFakingTau(UncertaintyScale scale, double eta,
                                                           GenLeptonMatch genLeptonMatch,
                                                           TauIdDiscriminator tauVSeDiscriminator, int decayMode) const
{
    std::vector<int> dms = {0, 1};
    std::vector<std::string> regions = {"barrel", "endcap"};
    std::map<std::pair<int, std::string>, double> fes; //mettere invece di string fare booleano
    std::map<std::pair<int, std::string>, std::pair<double, double>> fes_error;

    Int_t i = 0;
    bool is_barrel = (i == 0 || i == 1 ) ? true : false;
    for(const auto& region: regions){
        for(const auto& dm: dms){
            double y = hist_ele_faking_tau->GetY()[i];
            double y_error_up = hist_ele_faking_tau->GetErrorYhigh(i);
            double y_error_low = hist_ele_faking_tau->GetErrorYlow(i);
            fes[std::make_pair(dm, region)] = y;
            fes_error[std::make_pair(dm, region)] = std::make_pair(y_error_low, y_error_up);
            ++i;
        }
    }

    if((genLeptonMatch == GenLeptonMatch::Electron ||genLeptonMatch == GenLeptonMatch::TauElectron) &&
        tauVSeDiscriminator == TauIdDiscriminator::byDeepTau2017v2p1VSe &&
        (std::find(dms.begin(), dms.end(), decayMode) != dms.end())){

        std::string reg = "";
        if(std::abs(eta) < 1.5)
            reg = "barrel";
        else
            reg = "endcap";

        auto errors = fes_error[std::make_pair(decayMode, reg)];
        auto error_fes  = static_cast<int>(scale) > 0 ? errors.second : errors.first;
        auto fes_value = fes[std::make_pair(decayMode, reg)];

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
