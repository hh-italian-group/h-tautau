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
//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {
using TauIdDiscriminator = analysis::TauIdDiscriminator;
class TauESUncertainties{
public:
    static double GetCorrectionFactor(analysis::Period period, int decayMode, GenLeptonMatch genLeptonMatch,
                                      UncertaintySource unc_source, UncertaintyScale scale, double pt,
                                      TauIdDiscriminator tauVSjetDiscriminator, TauIdDiscriminator tauVSeDiscriminator,
                                      double eta);

    static double GetCorrectionFactorTrueTau(analysis::Period period, int decayMode, UncertaintyScale current_scale,
                                             double pt, TauIdDiscriminator tauVSjetDiscriminator);

    static double GetCorrectionFactorTrueMuon(analysis::Period period, int decayMode);

    static double GetCorrectionFactorEleFakingTau(analysis::Period period, UncertaintyScale scale, double eta,
                                                  TauIdDiscriminator tauVSeDiscriminator, int decayMode);
};
} // namespace analysis
