/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {
namespace sv_fit_ana {

struct FitResults {
    bool has_valid_momentum;
    LorentzVectorM momentum;
    LorentzVectorM momentum_error;
    double transverseMass;
    double transverseMass_error;

    FitResults() : has_valid_momentum(false), transverseMass(std::numeric_limits<double>::lowest()),
                   transverseMass_error(std::numeric_limits<double>::lowest()) {}
};

inline classic_svFit::MeasuredTauLepton CreateMeasuredLepton(const ntuple::TupleLepton& lepton);

class FitProducer {
public:
    explicit FitProducer(int _verbosity = 0)
        : verbosity(_verbosity)
    {
        TH1::AddDirectory(false);
    }

    FitResults Fit(const LeptonCandidate<ntuple::TupleLepton>& first_daughter,
                   const LeptonCandidate<ntuple::TupleLepton>& second_daughter,
                   const MissingET<ntuple::TupleMet>& met) const;

private:
    FitResults RunAlgorithm(const std::vector<classic_svFit::MeasuredTauLepton>& measured_leptons,
                            const LorentzVector& met_momentum, const SquareMatrix<2>& met_cov) const;

private:
    int verbosity;
};

} // namespace sv_fit_ana
} // namespace analysis
