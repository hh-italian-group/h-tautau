/*! Definition of wrapper to estimate jet energy uncertainties.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <TLorentzVector.h>

#include "FWCore/Utilities/src/EDMException.cc"
#include "FWCore/Utilities/src/Exception.cc"
#include "FWCore/Utilities/src/typelookup.cc"

#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc"

#include "TreeProduction/interface/Jet.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

class JetEnergyUncertaintyCorrector {
public:
    JetEnergyUncertaintyCorrector(const std::string& input_file_name, const std::string& section_name)
        : parameters(input_file_name, section_name), jetCorrector(parameters) {}

    void ApplyCorrection(ntuple::Jet& jet, bool scale_up)
    {
        try {
            jetCorrector.setJetPt(jet.pt);
            jetCorrector.setJetEta(jet.eta);
            const double uncertainty = jetCorrector.getUncertainty(scale_up);
            const double sign = scale_up ? +1 : -1;
            const double sf = 1.0 + sign * uncertainty;
            const TLorentzVector original_momentum = MakeLorentzVectorPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass);
            const TLorentzVector corrected_momentum = original_momentum * sf;
            jet.pt = corrected_momentum.Pt();
            jet.eta = corrected_momentum.Eta();
            jet.phi = corrected_momentum.Phi();
            jet.mass = corrected_momentum.M();
        } catch(cms::Exception&) {
            std::cerr << "WARNING: Jet uncertainty is not calculated for jet with pt = " << jet.pt
                      << ", eta = " << jet.eta << std::endl;
        }
    }

    void ApplyCorrection(ntuple::JetVector& jets, bool scale_up)
    {
        for(ntuple::Jet& jet : jets)
            ApplyCorrection(jet, scale_up);
    }

private:
    void CmsswWarningsWorkaround() const
    {
        float fx[0], fy[0];
        quadraticInterpolation(0, fx, fy);
    }

private:
    JetCorrectorParameters parameters;
    JetCorrectionUncertainty jetCorrector;
};

} // namespace analysis
