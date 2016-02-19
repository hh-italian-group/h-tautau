/*! Definition of wrappers for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "SVfit/source/generalAuxFunctions.cc"
#include "SVfit/source/LikelihoodFunctions.cc"
#include "SVfit/source/MarkovChainIntegrator.cc"
#include "SVfit/source/NSVfitStandaloneAlgorithm.cc"
#include "SVfit/source/NSVfitStandaloneLikelihood.cc"
#include "SVfit/source/svFitAuxFunctions.cc"

#include "TreeProduction/interface/MET.h"
#include "Candidate.h"

namespace analysis {

namespace sv_fit {

enum class FitAlgorithm { Vegas, MarkovChain };
std::ostream& operator<< (std::ostream& s, FitAlgorithm algo)
{
    static const std::map<FitAlgorithm, std::string> algo_names = { { FitAlgorithm::Vegas, "VEGAS" },
                                                                    { FitAlgorithm::MarkovChain, "Markov chain" } };
    s << algo_names.at(algo);
    return s;
}

struct FitResults {
    static constexpr double default_value = std::numeric_limits<double>::lowest();

    bool has_valid_mass;
    double mass;

    bool has_valid_momentum;
    TLorentzVector momentum;

    FitResults() : has_valid_mass(false), mass(default_value), has_valid_momentum(false) {}
};

struct CombinedFitResults {
    FitResults fit_vegas;
    FitResults fit_mc;
};

inline FitResults Fit(FitAlgorithm fitAlgorithm, const CandidatePtrVector& higgsLegs, const ntuple::MET& met)
{
    static const bool debug = false;

    static const std::map<Candidate::Type, NSVfitStandalone::kDecayType> decayTypeMap = {
        { Candidate::Type::Electron, NSVfitStandalone::kLepDecay },
        { Candidate::Type::Muon, NSVfitStandalone::kLepDecay },
        { Candidate::Type::Tau, NSVfitStandalone::kHadDecay }
    };

    if (higgsLegs.size() != 2)
        throw exception("Invalid number of legs to perform svFit");

    const Candidate higgsCandidate(Candidate::Type::Higgs, higgsLegs.at(0), higgsLegs.at(1));
    const NSVfitStandalone::Vector measuredMET(met.pt * std::cos(met.phi), met.pt * std::sin(met.phi), 0.0);
    const TMatrixD covMET = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);

    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    for (const auto& leg : higgsLegs){
        if(!decayTypeMap.count(leg->GetType()))
            throw exception("leg is not compatible with leptonic or hadronic tau decay");

        const NSVfitStandalone::kDecayType decayType = decayTypeMap.at(leg->GetType());
        const NSVfitStandalone::LorentzVector lepton(leg->GetMomentum().Px(), leg->GetMomentum().Py(),
                                                     leg->GetMomentum().Pz(), leg->GetMomentum().E());

        measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(decayType, lepton));
    }

    NSVfitStandalone::NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMET, 0);
    algo.addLogM(false);

    if(fitAlgorithm == FitAlgorithm::Vegas)
        algo.integrateVEGAS();
    else if(fitAlgorithm == FitAlgorithm::MarkovChain)
        algo.integrateMarkovChain();
    else
        throw exception("Unsupported algorithm.");

    FitResults result;

    if(algo.isValidSolution()) {
        result.mass = algo.mass();
        result.has_valid_mass = true;
        if(fitAlgorithm == FitAlgorithm::MarkovChain) {
            result.momentum.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.mass());
            result.has_valid_momentum = true;
        }
    } else
        std::cerr << "Can't fit with " << fitAlgorithm << std::endl;

    if(debug) {
        std::cout << std::fixed << std::setprecision(4)
                  << "\nOriginal mass = " << higgsCandidate.GetMomentum().M()
                  << "\nOriginal momentum = " << higgsCandidate.GetMomentum()
                  << "\nFirst daughter momentum = " << higgsCandidate.GetDaughters().at(0)->GetMomentum()
                  << "\nSecond daughter momentum = " << higgsCandidate.GetDaughters().at(1)->GetMomentum()
                  << "\nMET momentum = (" << met.pt << ", " << met.phi << ")"
                  << "\nMET covariance: " << covMET
                  << "\nSVfit algorithm = " << fitAlgorithm;
        if(result.has_valid_mass)
            std::cout << "\nSVfit mass = " << result.mass;
        if(result.has_valid_momentum)
            std::cout << "\nSVfit momentum = " << result.momentum;
        std::cout << std::endl;
    }

    return result;
}

CombinedFitResults CombinedFit(const CandidatePtrVector& higgsLegs, const ntuple::MET& met, bool fitWithVegas,
                               bool fitWithMarkovChain)
{
    CombinedFitResults result;
    if(fitWithVegas)
        result.fit_vegas = Fit(FitAlgorithm::Vegas, higgsLegs, met);
    if(fitWithMarkovChain)
        result.fit_mc = Fit(FitAlgorithm::MarkovChain, higgsLegs, met);
    return result;
}

} // namespace sv_fit

} // namespace analysis

