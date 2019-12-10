/*! Tau-related uncertainties.
If not specified otherwise, all definitions are taken from the TauID for 13 TeV TWiki:
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include "h-tautau/Analysis/include/TauUncertainties.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_energy_scale
namespace analysis {

    double GetCorrectionFactor(analysis::Period period, int decayMode, UncertaintyScale scale, double pt)
    {
        //put no shift and 2% of unc for missing DM
        //values taken from: https://indico.cern.ch/event/864131/contributions/3644021/attachments/1946837/3230164/Izaak_TauPOG_TauES_20191118.pdf
        static const std::map<analysis::Period, std::map<int, PhysicalValue>> tau_correction_factor = {
          { analysis::Period::Run2016, { {0, PhysicalValue(-0.1,0.7)},
                                         {1, PhysicalValue(-0.1,0.3)},
                                         {5, PhysicalValue(0.0,2.0)},
                                         {6, PhysicalValue(0.0,2.0)},
                                         {10, PhysicalValue(0.0,0.4)},
                                         {11, PhysicalValue(2.6,0.6)}    }},
          { analysis::Period::Run2017, { {0, PhysicalValue(-0.7,0.7)},
                                         {1, PhysicalValue(-1.1,0.3)},
                                         {5, PhysicalValue(0.0,2.0)},
                                         {6, PhysicalValue(0.0,2.0)},
                                         {10, PhysicalValue(0.5,0.5)},
                                         {11, PhysicalValue(1.7,0.6)}} },
         { analysis::Period::Run2018, { {0, PhysicalValue(-1.6,0.8)},
                                        {1, PhysicalValue(0.8,0.3)},
                                        {5, PhysicalValue(0.0,2.0)},
                                        {6, PhysicalValue(0.0,2.0)},
                                        {10, PhysicalValue(-0.9,0.4)},
                                        {11, PhysicalValue(1.3,1.0)}} }
        };

        double correction_factor;
        if(pt > 400)
            correction_factor = 1 + static_cast<int>(scale) * 0.03;
        else{
            if(!tau_correction_factor.count(period))
                throw exception("Period not found in tau correction map.");
            if(!tau_correction_factor.at(period).count(decayMode))
                throw exception("Decay mode not found in tau correction map.");
            PhysicalValue tau_correction = tau_correction_factor.at(period).at(decayMode);
            correction_factor = 1 + ((tau_correction.GetValue() + static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
            //double uncertainty = 1 + ((static_cast<int>(scale) * tau_correction.GetStatisticalError())/100);
        }
        return correction_factor;
    }

} // namespace analysis
