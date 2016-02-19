/*! Definition of wrapper for RecoilCorrector code to apply recoil correction.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <TLorentzVector.h>

#include "RecoilCorrector_v7/RecoilCorrector.hh"
#include "TreeProduction/interface/MET.h"

namespace analysis {

class RecoilCorrectionProducer{
public:
    RecoilCorrectionProducer(const std::string& fileCorrectTo, const std::string& fileZmmData,
                                 const std::string& fileZmmMC)
        : corrector(fileCorrectTo)
    {
        corrector.addDataFile(fileZmmData);
        corrector.addMCFile(fileZmmMC);
    }

    ntuple::MET ApplyCorrection(const ntuple::MET& originalMET, const TLorentzVector& resonantMomentum,
                                const TLorentzVector& resonantMomentumMC, size_t njets)
    {
        static const bool debug = false;
        //from Riccardo
    //    static const std::string fileCorrectTo = "RecoilCorrector_v7/recoilfits/recoilfit_htt53X_20pv_njet.root";
    //    static const std::string fileZmmData = "RecoilCorrector_v7/recoilfits/recoilfit_datamm53XRR_2012_njet.root";
    //    static const std::string fileZmmMC = "RecoilCorrector_v7/recoilfits/recoilfit_zmm53XRR_2012_njet.root";
        double iU1, iU2;
        const double iFluc = 0, iScale = 0;
        double met = originalMET.pt;
        double metphi = originalMET.phi;
    //    RecoilCorrector corrector(fileCorrectTo);
    //    corrector.addDataFile(fileZmmData);
    //    corrector.addMCFile(fileZmmMC);

        corrector.CorrectType1(met, metphi, resonantMomentumMC.Pt(), resonantMomentumMC.Phi(),
                               resonantMomentum.Pt(), resonantMomentum.Phi(), iU1, iU2, iFluc, iScale, njets);

        ntuple::MET correctedMET(originalMET);
        correctedMET.pt = met;
        correctedMET.phi = metphi;

        if(debug) {
            std::cout << "MET: " << originalMET.pt << ", MET phi: " << originalMET.phi
                      << "\nPt GEN: " << resonantMomentumMC.Pt() << ", Phi GEN: " << resonantMomentumMC.Phi()
                      << "\nPt RECO: " << resonantMomentum.Pt() << ", Phi RECO: " << resonantMomentum.Phi()
                      << "\njets: " << njets << ", Fluc: " << iFluc << ", Scale: " << iScale
                      << "\nCorrected MET: " << correctedMET.pt << ", Corrected MET phi: " << correctedMET.phi
                      << "\nU1: " << iU1 << ", U2: " << iU2 << std::endl;
        }

        return correctedMET;
    }

private:
    RecoilCorrector corrector;
};

} // analysis
