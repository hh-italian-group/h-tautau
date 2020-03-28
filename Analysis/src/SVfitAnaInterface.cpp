/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../include/SVfitAnaInterface.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

namespace analysis {

namespace sv_fit_ana {

classic_svFit::MeasuredTauLepton CreateMeasuredLepton(const LeptonCandidate<ntuple::TupleLepton>& lepton)
{
    const auto& momentum = lepton.GetMomentum();
    double preciseVisMass = momentum.mass();
    classic_svFit::MeasuredTauLepton::kDecayType decay_type;
    int decay_mode = -1;
    if(lepton->leg_type() == analysis::LegType::e) {
        static const double minVisMass = classic_svFit::electronMass, maxVisMass = minVisMass;
        preciseVisMass = std::clamp(preciseVisMass, minVisMass, maxVisMass);
        decay_type = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    } else if(lepton->leg_type() == analysis::LegType::mu) {
        decay_type = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    } else if(lepton->leg_type() == analysis::LegType::tau){
        const double minVisMass = lepton->decayMode() == 0 ? classic_svFit::chargedPionMass : 0.3;
        const double maxVisMass = lepton->decayMode() == 0 ? classic_svFit::chargedPionMass : 1.5;
        preciseVisMass = std::clamp(preciseVisMass, minVisMass, maxVisMass);
        decay_type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        decay_mode = lepton->decayMode();
    } else {
        throw exception("Leg Type not supported for SVFitAnaInterface.");
    }

    return classic_svFit::MeasuredTauLepton(decay_type, momentum.Pt(), momentum.Eta(), momentum.Phi(), preciseVisMass,
                                            decay_mode);
}

FitResults FitProducer::Fit(const LeptonCandidate<ntuple::TupleLepton>& first_daughter,
                            const LeptonCandidate<ntuple::TupleLepton>& second_daughter,
                            const MissingET<ntuple::TupleMet>& met, int verbosity)
{
    static const auto init = []() { TH1::AddDirectory(false); return true; };
    static const bool initialized = init();
    (void) initialized;

    const std::vector<classic_svFit::MeasuredTauLepton> measured_leptons = {
        CreateMeasuredLepton(first_daughter),
        CreateMeasuredLepton(second_daughter)
    };

    const TMatrixD met_cov_t = ConvertMatrix(met.GetCovMatrix());
    ClassicSVfit algo(verbosity);
    algo.addLogM_fixed(false);
    algo.addLogM_dynamic(false);
    // algo.setDiTauMassConstraint(-1.0);
    if(verbosity > 0) {
        std::cout << "SVfit inputs:\n";
        for(size_t n = 0; n < measured_leptons.size(); ++n) {
            const auto& lep = measured_leptons.at(n);
            std::cout << std::fixed << std::setprecision(7);
            std::cout << "\tlep" << n << ", (pt, eta, phi, m) = (" << lep.pt() << ", " << lep.eta()
                      << ", " << lep.phi() << ", " << lep.mass() << "), type=" << lep.type()
                      << ", decayMode=" << lep.decayMode() << "\n";
        }
        std::cout << "\tMET (px, py) = (" << met.GetMomentum().Px() << ", " << met.GetMomentum().Py() << ")\n"
                  << "\tmet_cov: (00, 01, 10, 11) = (" << met_cov_t[0][0] << ", " << met_cov_t[0][1]
                  << ", " << met_cov_t[1][0] << ", " << met_cov_t[1][1] << ")\n";
    }
    algo.integrate(measured_leptons, met.GetMomentum().Px(), met.GetMomentum().Py(), met_cov_t);

    FitResults result;
    if(algo.isValidSolution()) {
        auto histoAdapter = dynamic_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter());
        result.momentum = LorentzVectorM(histoAdapter->getPt(), histoAdapter->getEta(), histoAdapter->getPhi(), histoAdapter->getMass());
        result.momentum_error = LorentzVectorM(histoAdapter->getPtErr(), histoAdapter->getEtaErr(), histoAdapter->getPhiErr(), histoAdapter->getMassErr());
        result.transverseMass = histoAdapter->getTransverseMass();
        result.transverseMass_error = histoAdapter->getTransverseMassErr();
        result.has_valid_momentum = true;
        if(verbosity > 0) {
            std::cout << "SVfit result: (pt, eta, phi, m) = (" << histoAdapter->getPt() << ", "
                      << histoAdapter->getEta() << ", " << histoAdapter->getPhi() << ", "
                      << histoAdapter->getMass() << ")\n";
        }
    }
    return result;
}

} // namespace sv_fit
} // namespace analysis
