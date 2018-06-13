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

    template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename Met>
    FitResults Fit(const LVector1& lepton1_p4, const LVector2& lepton2_p4,
                   const LVector3& jet1_p4, const LVector4& jet2_p4, const Met& met,
                   double resolution_1, double resolution_2) const
    {
        const auto& met_p4 = met.GetMomentum();
        const TMatrixD& met_cov = ConvertMatrix(met.GetCovMatrix());

        FitResults result;
        try {
            if(verbosity >= 100) {
                std::cout << std::setprecision(6);
                std::cout << "lepton1_p4 = " << lepton1_p4 << std::endl;
                std::cout << "lepton2_p4 = " << lepton2_p4 << std::endl;
                std::cout << "jet1_p4 = " << jet1_p4 << std::endl;
                std::cout << "jet2_4 = " << jet2_p4 << std::endl;
                std::cout << "met_p4 = " << met_p4 << std::endl;
                std::cout << "met_cov:" << met_cov << std::endl;
            }

            // add here the resolution
            HHKinFit2::HHKinFitMasterHeavyHiggs hh_kin_fit(ConvertVector(jet1_p4), ConvertVector(jet2_p4),
                                                           ConvertVector(lepton1_p4), ConvertVector(lepton2_p4),
                                                           TVector2(met_p4.Px(), met_p4.Py()), met_cov,
                                                           resolution_1, resolution_2);
            hh_kin_fit.verbosity = verbosity;
            hh_kin_fit.fit();

            result.convergence = hh_kin_fit.getConvergence();
            if(result.HasValidMass()) {
                result.mass = hh_kin_fit.getMH();
                result.chi2 = hh_kin_fit.getChi2();
                result.probability = hh_kin_fit.getFitProb();
            }

            if(verbosity >= 100) {
                std::cout << "Convergence = " << result.convergence << std::endl;
                if(result.HasValidMass()) {
                    std::cout << std::setprecision(6);
                    std::cout << "Mass = " << result.mass << std::endl;
                    std::cout << "chi2 = " << result.chi2 << std::endl;
                    std::cout << "probability = " << result.probability << std::endl;
                }
            }
        } catch(std::exception&) {}

        return result;
    }

private:
    int verbosity;
};

} // namespace kin_fit
} // namespace analysis
