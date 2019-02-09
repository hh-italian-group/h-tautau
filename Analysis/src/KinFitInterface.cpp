/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {

namespace kin_fit {

FitProducer::FitProducer(int _verbosity) : verbosity(_verbosity) {}

FitResults FitProducer::FitImpl(const TLorentzVector& lepton1_p4, const TLorentzVector& lepton2_p4,
               const TLorentzVector& jet1_p4, const TLorentzVector& jet2_p4, const TVector2& met,
               const TMatrixD& met_cov, double resolution_1, double resolution_2) const
{
    FitResults result;
    try {
        if(verbosity >= 100) {
            const auto print_p4 = [](const std::string& name, const TLorentzVector& p4) {
                std::cout << name << " (pt, eta, phi, mass) = (" << p4.Pt() << ", " << p4.Eta() << ", "
                          << p4.Phi() << ", " << p4.M() << ")" << std::endl;
            };
            std::cout << std::setprecision(6);
            print_p4("lepton1", lepton1_p4);
            print_p4("lepton2", lepton2_p4);
            print_p4("jet1", jet1_p4);
            print_p4("jet2", jet2_p4);
            std::cout << "met (pt, phi) = (" << met.Mod() << ", " << met.Phi() << ")" << std::endl;
            std::cout << "met_cov:" << met_cov << std::endl;
        }

        // add here the resolution
        HHKinFit2::HHKinFitMasterHeavyHiggs hh_kin_fit(jet1_p4, jet2_p4, lepton1_p4, lepton2_p4, met, met_cov,
                                                       resolution_1, resolution_2);
        hh_kin_fit.verbosity = verbosity;
        hh_kin_fit.fit();

        result.convergence = hh_kin_fit.getConvergence();
        if(result.HasValidMass()) {
            result.mass = hh_kin_fit.getMH();
            result.chi2 = hh_kin_fit.getChi2();
            result.probability = hh_kin_fit.getFitProb();
        }

        if(verbosity >= 100) {
            std::cout << "Convergence = " << result.convergence << std::endl;
            if(result.HasValidMass()) {
                std::cout << std::setprecision(6);
                std::cout << "Mass = " << result.mass << std::endl;
                std::cout << "chi2 = " << result.chi2 << std::endl;
                std::cout << "probability = " << result.probability << std::endl;
            }
        }
    } catch(std::exception&) {}

    return result;
}

} // namespace kin_fit
} // namespace analysis
