/*! Definition of wrappers for HHKinFit code.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "HHKinFit/src/HHDiJetKinFit.cpp"
#include "HHKinFit/src/HHDiJetKinFitMaster.cpp"
#include "HHKinFit/src/HHEventRecord.cpp"
#include "HHKinFit/src/HHKinFit.cpp"
#include "HHKinFit/src/HHKinFitMaster.cpp"
#include "HHKinFit/src/HHParticle.cpp"
#include "HHKinFit/src/HHParticleList.cpp"
//#include "HHKinFit/src/HHTauInputKinFit.cpp"
//#include "HHKinFit/src/HHTauInputKinFitMaster.cpp"
#include "HHKinFit/src/HHV4Vector.cpp"
#include "HHKinFit/src/PSMath.cpp"
#include "HHKinFit/src/PSTools.cpp"

#include "TreeProduction/interface/MET.h"
#include "AnalysisBase/include/Candidate.h"
#include "AnalysisBase/include/RootExt.h"

namespace analysis {

namespace kinematic_fit {

namespace four_body {
struct FitInput {
    std::vector<TLorentzVector> bjet_momentums;
    std::vector<TLorentzVector> tau_momentums;
    TLorentzVector mvaMET;
    TMatrixD metCov;

    FitInput(const TLorentzVector& bjet1, const TLorentzVector& bjet2,
                     const TLorentzVector& tau1, const TLorentzVector& tau2,
                     const TLorentzVector& _mvaMET, const TMatrixD& _metCov)
        : mvaMET(_mvaMET), metCov(_metCov)
    {
        bjet_momentums.push_back(bjet1);
        bjet_momentums.push_back(bjet2);
        tau_momentums.push_back(tau1);
        tau_momentums.push_back(tau2);
    }

    FitInput(const FitInput& input, double bjet_scale_factor, double tau_scale_factor)
        : mvaMET(input.mvaMET), metCov(input.metCov)
    {
        for(const TLorentzVector& momentum : input.bjet_momentums)
            bjet_momentums.push_back(momentum * bjet_scale_factor);

        for(const TLorentzVector& momentum : input.tau_momentums)
            tau_momentums.push_back(momentum * tau_scale_factor);
    }
};

struct FitResults {
    static constexpr double default_double_value = std::numeric_limits<double>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    bool has_valid_mass;
    double mass;
    int convergence;
    double chi2;
    double fit_probability;
    double pull_balance;
    double pull_balance_1;
    double pull_balance_2;

    FitResults()
        : has_valid_mass(false), mass(default_double_value), convergence(default_int_value), chi2(default_double_value),
          fit_probability(default_double_value), pull_balance(default_double_value),
          pull_balance_1(default_double_value), pull_balance_2(default_double_value)
    {}
};

inline FitResults Fit(const FitInput& input)
{
    static const bool debug = false;
    static const Int_t higgs_mass_hypotesis = 125;
//    static const Int_t convergence_cut = 0;
//    static const Double_t chi2_cut = 25;
//    static const Double_t pull_balance_cut = 0;

    const std::vector<Int_t> hypo_mh1 = { higgs_mass_hypotesis };
    const std::vector<Int_t> hypo_mh2 = { higgs_mass_hypotesis };

    if(debug) {
        std::cout << "Format: (Pt, eta, phi, E)\n";
        std::cout << "b1 momentum: " << input.bjet_momentums.at(0) << std::endl;
        std::cout << "b2 momentum: " << input.bjet_momentums.at(1) << std::endl;
        std::cout << "tau1 momentum: " << input.tau_momentums.at(0) << std::endl;
        std::cout << "tau2 momentum: " << input.tau_momentums.at(1) << std::endl;
        std::cout << "MET: " << input.mvaMET << std::endl;
        std::cout << "MET covariance:\n";
        input.metCov.Print();
        std::cout << "metDet = " << input.metCov.Determinant() << std::endl;
    }

    //intance of fitter master class
    TLorentzVector b1(input.bjet_momentums.at(0)), b2(input.bjet_momentums.at(1)),
            tau1(input.tau_momentums.at(0)), tau2(input.tau_momentums.at(1)),
            mvaMET(input.mvaMET);
    HHKinFitMaster kinFit(&b1, &b2, &tau1, &tau2);
    kinFit.setAdvancedBalance(&mvaMET, input.metCov);
    //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
    kinFit.addMh1Hypothesis(hypo_mh1);
    kinFit.addMh2Hypothesis(hypo_mh2);
    kinFit.doFullFit();

    if(debug) {
        const Double_t chi2_best = kinFit.getBestChi2FullFit();
        const Double_t mh_best = kinFit.getBestMHFullFit();
        const std::pair<Int_t, Int_t> bestHypo = kinFit.getBestHypoFullFit();

        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
        std::cout << "best mH=         " << mh_best << std::endl;
    }

    const std::pair<Int_t, Int_t> hypo(hypo_mh1.at(0), hypo_mh2.at(0));
    const Int_t convergence = kinFit.getConvergenceFullFit().at(hypo);
    const Double_t chi2 = kinFit.getChi2FullFit().at(hypo);
    const Double_t fitprob = kinFit.getFitProbFullFit().at(hypo);
    const Double_t pull_b1 = kinFit.getPullB1FullFit().at(hypo);
    const Double_t pull_b2 = kinFit.getPullB2FullFit().at(hypo);
    const Double_t pull_balance = kinFit.getPullBalanceFullFit().at(hypo);
    const Double_t mH = kinFit.getMHFullFit().at(hypo);

    if (debug) {
        std::cout << "fit convergence =  " << convergence << std::endl;
        std::cout << "fit chi2 =         " << chi2 << std::endl;
        std::cout << "fit fitprob =      " << fitprob << std::endl;
        std::cout << "fit pull_b1 =      " << pull_b1 << std::endl;
        std::cout << "fit pull_b2 =      " << pull_b2 << std::endl;
        std::cout << "fit pull_balance = " << pull_balance << std::endl;
        std::cout << "fit mH =           " << mH << std::endl;
    }

    FitResults result;
    result.convergence = convergence;
    result.chi2 = chi2;
    result.mass = mH;
    result.fit_probability = fitprob;
    result.pull_balance = pull_balance;
    result.pull_balance_1 = pull_b1;
    result.pull_balance_2 = pull_b2;
    result.has_valid_mass = convergence > 0;
    if (!result.has_valid_mass)
        result.mass = (b1+b2+tau1+tau2+mvaMET).M();
    if(!result.has_valid_mass && debug)
        std::cerr << "four body mass with kin Fit cannot be calculated" << std::endl;

    return result;
}

} // namespace four_body

namespace two_body {

struct FitInput {
    std::vector<TLorentzVector> bjet_momentums;

    FitInput(const TLorentzVector& bjet1, const TLorentzVector& bjet2)
    {
        bjet_momentums.push_back(bjet1);
        bjet_momentums.push_back(bjet2);
    }

    FitInput(const FitInput& input, double bjet_scale_factor)
    {
        for(const TLorentzVector& momentum : input.bjet_momentums)
            bjet_momentums.push_back(momentum * bjet_scale_factor);
    }
};

struct FitResults {
    static constexpr double default_double_value = std::numeric_limits<double>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    int convergence;
    double chi2;
    std::vector<TLorentzVector> bjet_momentums;

    FitResults() : convergence(default_int_value), chi2(default_double_value) {}
};


FitResults Fit(const FitInput& input)
{
    static const bool debug = false;
    static const Int_t higgs_mass_hypotesis = 125;

    const std::vector<Int_t> hypo_mh = { higgs_mass_hypotesis };

    if(debug) {
        std::cout << "Format: (Pt, eta, phi, E)\n";
        std::cout << "b1 momentum: " << input.bjet_momentums.at(0) << std::endl;
        std::cout << "b2 momentum: " << input.bjet_momentums.at(1) << std::endl;
    }

    //intance of fitter master class
    TLorentzVector b1(input.bjet_momentums.at(0)), b2(input.bjet_momentums.at(1));
    HHDiJetKinFitMaster kinFit(&b1, &b2);
    kinFit.addMhHypothesis(hypo_mh);
    kinFit.doFullFit();

    if(debug) {
        const Double_t chi2_best = kinFit.getBestChi2FullFit();
        const std::pair<Int_t, Int_t> bestHypo = kinFit.getBestHypoFullFit();

        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
    }

    const std::pair<Int_t, Int_t> hypo(hypo_mh.at(0), -1);
    const Int_t convergence = kinFit.getConvergenceFullFit().at(hypo);
    const Double_t chi2 = kinFit.getBestChi2FullFit();

    if (debug) {
        const Double_t fitprob = kinFit.getFitProbFullFit().at(hypo);
        const Double_t pull_b1 = kinFit.getPullB1FullFit().at(hypo);
        const Double_t pull_b2 = kinFit.getPullB2FullFit().at(hypo);

        std::cout << "fit convergence =  " << convergence << std::endl;
        std::cout << "fit chi2 =         " << chi2 << std::endl;
        std::cout << "fit fitprob =      " << fitprob << std::endl;
        std::cout << "fit pull_b1 =      " << pull_b1 << std::endl;
        std::cout << "fit pull_b2 =      " << pull_b2 << std::endl;
    }

    FitResults result;
    result.convergence = convergence;
    result.chi2 = chi2;
    result.bjet_momentums.push_back(kinFit.getFitJet1());
    result.bjet_momentums.push_back(kinFit.getFitJet2());

    return result;
}

} // namespace two_body

typedef std::pair<size_t, size_t> FitId;
typedef std::map<FitId, four_body::FitResults> FitResultsMap;

} // namespace kinematic_fit
} // namespace analysis
