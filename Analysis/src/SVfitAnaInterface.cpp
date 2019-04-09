/*! A wrapper for SVfit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../include/SVfitAnaInterface.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

namespace analysis {

namespace sv_fit_ana {

classic_svFit::MeasuredTauLepton CreateMeasuredLepton(const ntuple::TupleLepton& lepton)
{
    if(lepton.leg_type() == analysis::LegType::e){
        const auto& momentum = lepton.p4();
        // applying fix for electron mass
        static const double minVisMass = classic_svFit::electronMass, maxVisMass = minVisMass;
        double preciseVisMass = momentum.mass();
        if ( preciseVisMass < minVisMass ) preciseVisMass = minVisMass;
        if ( preciseVisMass > maxVisMass ) preciseVisMass = maxVisMass;
        return classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay,
                                                momentum.Pt(), momentum.Eta(), momentum.Phi(), preciseVisMass);
    }
    if(lepton.leg_type() == analysis::LegType::mu){
        const auto& momentum = lepton.p4();
        return classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay,
                                                  momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.M());
    }
    if(lepton.leg_type() == analysis::LegType::tau){
        const auto& momentum = lepton.p4();
        return classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,
                                                  momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.M(),
                                                  lepton.decayMode());
    }
    throw exception("Leg Type not supported for SVFitAnaInterface.");
}

FitProducer::FitProducer(int _verbosity)
    : verbosity(_verbosity)
{
    TH1::AddDirectory(false);
}

FitResults FitProducer::Fit(const LeptonCandidate<ntuple::TupleLepton>& first_daughter,
                            const LeptonCandidate<ntuple::TupleLepton>& second_daughter,
                            const MissingET<ntuple::TupleMet>& met) const
{
    std::vector<classic_svFit::MeasuredTauLepton> measured_leptons = {
        CreateMeasuredLepton(*first_daughter),
        CreateMeasuredLepton(*second_daughter)
    };

    const TMatrixD met_cov_t = ConvertMatrix(met.GetCovMatrix());
    ClassicSVfit algo(verbosity);
    algo.addLogM_fixed(false);
    algo.addLogM_dynamic(false);
    // algo.setDiTauMassConstraint(-1.0);
    algo.integrate(measured_leptons, met.GetMomentum().Px(), met.GetMomentum().Py(), met_cov_t);

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
