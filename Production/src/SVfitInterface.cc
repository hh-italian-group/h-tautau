/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/SVfitInterface.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

namespace analysis {

namespace sv_fit {

FitResults FitProducer::RunAlgorithm(const std::vector<svFitStandalone::MeasuredTauLepton>& measured_leptons,
                                     const LorentzVector& met_momentum, const SquareMatrix<2>& met_cov) const
{
    const TMatrixD met_cov_t = ConvertMatrix(met_cov);
    SVfitStandaloneAlgorithm algo(measured_leptons, met_momentum.Px(), met_momentum.Py(), met_cov_t, verbosity);
    algo.addLogM(false);
    algo.shiftVisPt(true, visPtResolutionFile.get());
    algo.integrateMarkovChain();

    FitResults result;
    if(algo.isValidSolution()) {
        result.momentum = LorentzVectorM(algo.pt(), algo.eta(), algo.phi(), algo.mass());
        result.transverseMass = algo.transverseMass();
        result.has_valid_momentum = true;
    }
    return result;
}

} // namespace sv_fit
} // namespace analysis
