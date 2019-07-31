/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include <utility>
#include <string>
#include <iostream>

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

    double GetCorrectionFactor(analysis::Period period, int decayMode, UncertaintyScale scale, double pt);

} // namespace analysis
