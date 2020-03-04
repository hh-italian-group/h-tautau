/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include <utility>
#include <string>
#include <iostream>
#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"


//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {
using TauIdDiscriminator = analysis::TauIdDiscriminator;
class TauESUncertainties{
public:
    TauESUncertainties(std::string file_low_pt, std::string file_high_pt,
                       DiscriminatorWP _ele_id_wp, std::string files_ele_faking_tau);
    double GetCorrectionFactor(int decayMode, GenLeptonMatch genLeptonMatch,
                               UncertaintySource unc_source, UncertaintyScale scale, double pt,
                               TauIdDiscriminator tauVSeDiscriminator, double eta);

    double GetCorrectionFactorTrueTau(double pt, int decayMode, UncertaintyScale scale,
                                      GenLeptonMatch genLeptonMatch = GenLeptonMatch::Tau);

    double GetCorrectionFactorMuonFakingTau(UncertaintyScale scale);

    double GetCorrectionFactorEleFakingTau(UncertaintyScale scale, double eta, GenLeptonMatch genLeptonMatch, TauIdDiscriminator tauVSeDiscriminator,
                                           int decayMode);

private:
    DiscriminatorWP ele_id_wp;
    std::shared_ptr<TFile> file_low;
    std::shared_ptr<TH1F> hist_tes_pt_low;
    std::shared_ptr<TFile> file_high;
    std::shared_ptr<TH1F> hist_tes_pt_high;
    std::shared_ptr<TFile> file_ele_faking_tau;
    std::shared_ptr<TH1F> hist_ele_faking_tau;
};
} // namespace analysis
