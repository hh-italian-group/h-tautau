/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

namespace kin_fit {

struct FitResults {
    double mass, chi2, probability;
    int convergence;
    bool HasValidMass() const { return convergence > 0; }

    FitResults() : convergence(std::numeric_limits<int>::lowest()) {}
    FitResults(double _mass, double _chi2, double _probability, int _convergence) :
                mass(_mass), chi2(_chi2), probability(_probability), convergence(_convergence) {}
};

class FitProducer {
public:
    explicit FitProducer(int _verbosity = 0);

    template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename Met>
    FitResults Fit(const LVector1& lepton1_p4, const LVector2& lepton2_p4,
                   const LVector3& jet1_p4, const LVector4& jet2_p4, const Met& met,
                   double resolution_1, double resolution_2) const
    {
        const auto& met_p4 = met.GetMomentum();
        const TMatrixD& met_cov = ConvertMatrix(met.GetCovMatrix());
        return FitImpl(ConvertVector(lepton1_p4), ConvertVector(lepton2_p4), ConvertVector(jet1_p4),
                       ConvertVector(jet2_p4), TVector2(met_p4.Px(), met_p4.Py()), met_cov,
                       resolution_1, resolution_2);
    }

private:
    FitResults FitImpl(const TLorentzVector& lepton1_p4, const TLorentzVector& lepton2_p4,
                       const TLorentzVector& jet1_p4, const TLorentzVector& jet2_p4, const TVector2& met,
                       const TMatrixD& met_cov, double resolution_1, double resolution_2) const;

private:
    int verbosity;
};

} // namespace kin_fit
} // namespace analysis
