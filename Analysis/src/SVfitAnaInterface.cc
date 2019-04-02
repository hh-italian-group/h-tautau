/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../include/SVfitAnaInterface.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

namespace analysis {

namespace sv_fit_ana {

FitResults FitProducer::RunAlgorithm(const std::vector<classic_svFit::MeasuredTauLepton>& measured_leptons,
                                     const LorentzVector& met_momentum, const SquareMatrix<2>& met_cov) const
{
    const TMatrixD met_cov_t = ConvertMatrix(met_cov);
    ClassicSVfit algo(verbosity);
    algo.addLogM_fixed(false);
    algo.addLogM_dynamic(false);
    // algo.setDiTauMassConstraint(-1.0);
    algo.integrate(measured_leptons, met_momentum.Px(), met_momentum.Py(), met_cov_t);

    FitResults result;
    if(algo.isValidSolution()) {
        auto histoAdapter = dynamic_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter());
        result.momentum = LorentzVectorM(histoAdapter->getPt(), histoAdapter->getEta(), histoAdapter->getPhi(), histoAdapter->getMass());
        result.momentum_error = LorentzVectorM(histoAdapter->getPtErr(), histoAdapter->getEtaErr(), histoAdapter->getPhiErr(), histoAdapter->getMassErr());
        result.transverseMass = histoAdapter->getTransverseMass();
        result.transverseMass_error = histoAdapter->getTransverseMassErr();
        result.has_valid_momentum = true;
    }
    return result;
}

} // namespace sv_fit
} // namespace analysis
