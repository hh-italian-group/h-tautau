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
    static double GetCorrectionFactor(analysis::Period period, int decayMode, GenLeptonMatch genLeptonMatch,
                                      UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                      TauIdDiscriminator tauVSeDiscriminator,
                                      double eta, std::string file_low_pt, std::string file_high_pt);

    static double GetCorrectionFactorTrueTau(double pt, int decayMode, std::string file_low_pt,
                                             std::string file_high_pt, UncertaintyScale scale,
                                             GenLeptonMatch genLeptonMatch = GenLeptonMatch::Tau,
                                             UncertaintySource unc_source = UncertaintySource::None);

    static double GetCorrectionFactorMuonFakingTau(analysis::Period period, int decayMode);

    static double GetCorrectionFactorEleFakingTau(analysis::Period period, UncertaintyScale scale, double eta,
                                                  TauIdDiscriminator tauVSeDiscriminator, int decayMode);
};
} // namespace analysis
