/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"

#include "h-tautau/Analysis/include/Candidate.h"
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {
namespace sv_fit {

struct FitResults {
    bool has_valid_momentum;
    LorentzVectorM momentum;
    double transverseMass;

    FitResults() : has_valid_momentum(false), transverseMass(std::numeric_limits<double>::lowest()) {}
};

namespace detail {
template<typename Lepton>
inline svFitStandalone::MeasuredTauLepton CreateMeasuredLepton(const Lepton& lepton);

template<>
inline svFitStandalone::MeasuredTauLepton CreateMeasuredLepton(
        const LeptonCandidate<pat::Electron, edm::Ptr<pat::Electron>>& lepton)
{
    const auto& momentum = lepton.GetMomentum();
    return svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay,
                                              momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.M());
}

template<>
inline svFitStandalone::MeasuredTauLepton CreateMeasuredLepton(const LeptonCandidate<pat::Muon>& lepton)
{
    const auto& momentum = lepton.GetMomentum();
    return svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay,
                                              momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.M());
}

template<>
inline svFitStandalone::MeasuredTauLepton CreateMeasuredLepton(const LeptonCandidate<pat::Tau>& lepton)
{
    const auto& momentum = lepton.GetMomentum();
    return svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay,
                                              momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.M(),
                                              lepton->decayMode());
}

} // namespace detail

class FitProducer {
public:
    explicit FitProducer(const std::string& visPtResolutionFileName, int _verbosity = 0)
        : visPtResolutionFile(root_ext::OpenRootFile(visPtResolutionFileName)), verbosity(_verbosity)
    {
        TH1::AddDirectory(false);
    }

    template<typename FirstLeg, typename SecondLeg, typename MetObject>
    FitResults Fit(const CompositCandidate<FirstLeg, SecondLeg>& higgs, const MissingET<MetObject>& met) const
    {
        std::vector<svFitStandalone::MeasuredTauLepton> measured_leptons = {
            detail::CreateMeasuredLepton(higgs.GetFirstDaughter()),
            detail::CreateMeasuredLepton(higgs.GetSecondDaughter())
        };

        return RunAlgorithm(measured_leptons, met.GetMomentum(), met.GetCovMatrix());
    }

private:
    FitResults RunAlgorithm(const std::vector<svFitStandalone::MeasuredTauLepton>& measured_leptons,
                            const LorentzVector& met_momentum, const SquareMatrix<2>& met_cov) const;

private:
    std::shared_ptr<TFile> visPtResolutionFile;
    int verbosity;
};

} // namespace sv_fit
} // namespace analysis
