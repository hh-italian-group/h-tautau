/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/SVfitInterface.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

namespace analysis {

namespace sv_fit {

FitResults FitProducer::RunAlgorithm(const std::vector<classic_svFit::MeasuredTauLepton>& measured_leptons,
                                     const LorentzVector& met_momentum, const SquareMatrix<2>& met_cov) const
{
    const TMatrixD met_cov_t = ConvertMatrix(met_cov);
    ClassicSVfit algo(verbosity);
    algo.addLogM_fixed(false);
    // algo.setDiTauMassConstraint(-1.0);
    algo.integrate(measured_leptons, met_momentum.Px(), met_momentum.Py(), met_cov_t);

    FitResults result;
    if(algo.isValidSolution()) {
        auto* histoAdapter = dynamic_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter());
        double pt = histoAdapter->getPt();
        double eta = histoAdapter->getEta();
        double phi = histoAdapter->getPhi();
        double mass = histoAdapter->getMass();
        result.momentum = LorentzVectorM(pt, eta, phi, mass);
        double transverseMass = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass();
        result.transverseMass = transverseMass;
        result.has_valid_momentum = true;
    }
    return result;
}

} // namespace sv_fit
} // namespace analysis
