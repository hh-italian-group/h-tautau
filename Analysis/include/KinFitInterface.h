/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

namespace kin_fit {

struct FitResults {
    double mass, chi2, probability;
    int convergence;
    bool HasValidMass() const { return convergence > 0; }

    FitResults() : convergence(std::numeric_limits<int>::lowest()) {}
};

class FitProducer {
public:
    explicit FitProducer(int _verbosity = 0)
        : verbosity(_verbosity)
    {
    }

    template<typename Met>
    FitResults Fit(const std::vector<LorentzVector>& lepton_momentums,
                   const std::vector<LorentzVector>& jet_momentums, const Met& met)
    {
        if(lepton_momentums.size() != 2)
            throw exception("Invalid number of leptons to fit.");
        if(jet_momentums.size() != 2)
            throw exception("Invalid number of jets to fit.");

        const auto& met_momentum = met.GetMomentum();
        const TMatrixD met_cov = ConvertMatrix(met.GetCovMatrix());

        HHKinFit2::HHKinFitMasterHeavyHiggs hh_kin_fit(ConvertVector(lepton_momentums.at(0)),
                                                       ConvertVector(lepton_momentums.at(1)),
                                                       ConvertVector(jet_momentums.at(0)),
                                                       ConvertVector(jet_momentums.at(1)),
                                                       TVector2(met_momentum.Px(), met_momentum.Py()), met_cov);
        hh_kin_fit.verbosity = verbosity;
        hh_kin_fit.fit();
        FitResults result;
        result.convergence = hh_kin_fit.getConvergence();
        if(result.HasValidMass()) {
            result.mass = hh_kin_fit.getMH();
            result.chi2 = hh_kin_fit.getChi2();
            result.probability = hh_kin_fit.getFitProb();
        }
        return result;
    }

private:
    int verbosity;
};

} // namespace kin_fit
} // namespace analysis
