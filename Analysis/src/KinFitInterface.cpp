/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

namespace kin_fit {

FitResults FitProducer::FitImpl(const TLorentzVector& lepton1_p4, const TLorentzVector& lepton2_p4,
               const TLorentzVector& jet1_p4, const TLorentzVector& jet2_p4, const TVector2& met,
               const TMatrixD& met_cov, double resolution_1, double resolution_2, int verbosity)
{
    FitResults result;
    try {
        if(verbosity > 0) {
            std::cout << std::setprecision(6)
                      << "lepton1 " << analysis::LorentzVectorToString(lepton1_p4, analysis::LVectorRepr::PtEtaPhiME)
                      << "\nlepton2 " << analysis::LorentzVectorToString(lepton2_p4, analysis::LVectorRepr::PtEtaPhiME)
                      << "\njet1 "  << analysis::LorentzVectorToString(jet1_p4, analysis::LVectorRepr::PtEtaPhiME)
                      << ", resolution=" << resolution_1
                      << "\njet2 "  << analysis::LorentzVectorToString(jet2_p4, analysis::LVectorRepr::PtEtaPhiME)
                      << ", resolution=" << resolution_2
                      << "\nmet (px, py, pt, phi) = (" << met.Px() << ", " << met.Py() << ", "
                      << met.Mod() << ", " << met.Phi() << ")"
                      << "\nmet_cov: (00, 01, 10, 11) = (" << met_cov[0][0] << ", " << met_cov[0][1]
                      << ", " << met_cov[1][0] << ", " << met_cov[1][1] << ")\n";
        }

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

        if(verbosity > 0) {
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
