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

    template<typename LorentzVector1, typename LorentzVector2, typename Met>
    FitResults Fit(const std::vector<LorentzVector1>& lepton_momentums,
                   const std::vector<LorentzVector2>& jet_momentums, const Met& met)
    {
        if(lepton_momentums.size() != 2)
            throw exception("Invalid number of leptons to fit.");
        if(jet_momentums.size() != 2)
            throw exception("Invalid number of jets to fit.");

        const auto& met_momentum = met.GetMomentum();
        const TMatrixD met_cov = ConvertMatrix(met.GetCovMatrix());

        FitResults result;
        try {
            if(verbosity > 100) {
                std::cout << "l1_p4 = " << lepton_momentums.at(0) << "l2_p4 = " << lepton_momentums.at(1) << std::endl;
                std::cout << "b1_p4 = " << jet_momentums.at(0) << "b2_4 = " << jet_momentums.at(1) << std::endl;
                std::cout << "met_p4 = " << met_momentum << std::endl;
            }

            HHKinFit2::HHKinFitMasterHeavyHiggs hh_kin_fit(ConvertVector(jet_momentums.at(0)),
                                                           ConvertVector(jet_momentums.at(1)),
                                                           ConvertVector(lepton_momentums.at(0)),
                                                           ConvertVector(lepton_momentums.at(1)),
                                                           TVector2(met_momentum.Px(), met_momentum.Py()), met_cov);
            hh_kin_fit.verbosity = verbosity;
            hh_kin_fit.fit();

            result.convergence = hh_kin_fit.getConvergence();
            if(result.HasValidMass()) {
                result.mass = hh_kin_fit.getMH();
                result.chi2 = hh_kin_fit.getChi2();
                result.probability = hh_kin_fit.getFitProb();
            }
        } catch(std::exception&) {}

        return result;
    }

private:
    int verbosity;
};

} // namespace kin_fit
} // namespace analysis
