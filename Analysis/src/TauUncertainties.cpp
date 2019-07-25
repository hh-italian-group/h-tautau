/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

    std::pair<double,double> GetCorrectionFactor(analysis::Period period, int decayMode, UncertaintyScale scale)
    {
        static const std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor = {
          { analysis::Period::Run2016, { {0, PhysicalValue(-0.6,1.0)},
                                         {1, PhysicalValue(-0.5,0.9)},
                                         {10, PhysicalValue(0.0,1.1)} }},
          { analysis::Period::Run2017, { {0, PhysicalValue(0.7,0.8)},
                                         {1, PhysicalValue(-0.2,0.8)},
                                         {10, PhysicalValue(0.1,0.9)} } },
         { analysis::Period::Run2018, { {0, PhysicalValue(-1.3,1.1)},
                                        {1, PhysicalValue(-0.5,0.9)},
                                        {10, PhysicalValue(-1.2,0.8)} } }
        };

        std::pair<double,double> correction_factor;
        PhysicalValue tau_correction = tau_correction_factor.at(period).at(decayMode);
        double result = 1 + ((tau_correction.GetValue() + static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
        double uncertainty = 1 + ((static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
        correction_factor.first = result;
        correction_factor.second = uncertainty;
        return correction_factor;
    }

} // namespace analysis
