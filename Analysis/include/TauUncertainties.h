/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <TH1.h>
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {
class TauESUncertainties{
public:
    TauESUncertainties(const std::string& file_low_pt, const std::string& file_high_pt,
                       const std::string& file_ele_faking_tau);

    double GetCorrectionFactor(int decayMode, GenLeptonMatch genLeptonMatch,
                               UncertaintySource unc_source, UncertaintyScale scale, double pt, double eta) const;
    double GetCorrectionFactorTrueTau(double pt, int decayMode, UncertaintyScale scale) const;
    double GetCorrectionFactorEleFakingTau(double eta, int decayMode, UncertaintyScale scale) const;
    double GetCorrectionFactorMuonFakingTau(UncertaintyScale scale) const;

private:
    static std::map<int, StVariable> LoadTauCorrections(const std::string& file_name);
    static std::map<std::pair<int, bool>, StVariable> LoadElectronCorrections(const std::string& file_name);
    static bool ApplyUncertaintyScale(int decayMode, GenLeptonMatch genLeptonMatch, UncertaintySource unc_source);

private:
    std::map<int, StVariable> tes_low_pt, tes_high_pt;
    std::map<std::pair<int, bool>, StVariable> fes;
};
} // namespace analysis
