/*! Apply jet uncertainties to the event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <string>
#include <vector>
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "JetCorrectorParameters.h"   // CondFormats/JetMETObjects/interface
#include "JetCorrectionUncertainty.h" // CondFormats/JetMETObjects/interface

namespace jec {

using JetCandidate = analysis::Candidate<ntuple::TupleJet>;
using JetCollection = std::vector<JetCandidate>;
using analysis::UncertaintySource;
using analysis::UncertaintyScale;

class JECUncertaintiesWrapper
{
public:
    static const std::set<UncertaintySource>& JetFullUncertainties();
    static const std::set<UncertaintySource>& JetReducedUncertainties();

    JECUncertaintiesWrapper(const std::string& uncertainties_source, bool is_full, analysis::Period& period);

    const std::string ReturnJecName(UncertaintySource unc_source, bool is_full, analysis::Period& period);
    static bool IsJetUncertainties(UncertaintySource unc_source);

    template<typename JetCollection, typename LorentzVector1 = analysis::LorentzVector,
             typename LorentzVector2 = analysis::LorentzVector>
    JetCollection ApplyShift(const JetCollection& jet_candidates,
        analysis::UncertaintySource uncertainty_source,
        analysis::UncertaintyScale scale,
        const std::vector<LorentzVector1>* other_jets_p4 = nullptr,
        LorentzVector2* met = nullptr) const
    {
         static const std::map<UncertaintyScale, bool> scales = {
             { UncertaintyScale::Up, true }, { UncertaintyScale::Down, false }
         };

         static const std::map<UncertaintyScale, int> scales_variation = {
             { UncertaintyScale::Up, +1 }, { UncertaintyScale::Down, -1 }
         };

        JetCollection corrected_jets;
         if(!uncertainty_map.count(uncertainty_source))
             throw analysis::exception("Jet Uncertainty source % not found.") % uncertainty_source;
         if(scale == analysis::UncertaintyScale::Central)
             throw analysis::exception("Uncertainty scale Central.");
         auto unc = uncertainty_map.at(uncertainty_source);
         double shifted_met_px = 0;
         double shifted_met_py = 0;

         for (const auto& jet : jet_candidates){
             unc->setJetPt(static_cast<float>(jet.GetMomentum().pt()));
             unc->setJetEta(static_cast<float>(jet.GetMomentum().eta()));
             const auto unc_result = unc->getUncertainty(scales.at(scale));
             const double unc_var = unc_result ? *unc_result : 0;
             const int sign = scales_variation.at(scale);
             const auto sf = static_cast<typename LorentzVector1::Scalar>(1.0 + (sign * unc_var));
             const auto shiftedMomentum = jet.GetMomentum() * sf;
             auto corr_jet(jet);
             corr_jet.SetMomentum(shiftedMomentum);
             corrected_jets.push_back(corr_jet);
             shifted_met_px += jet.GetMomentum().px() - corr_jet.GetMomentum().px();
             shifted_met_py += jet.GetMomentum().py() - corr_jet.GetMomentum().py();
         }

         if(met){
            if(other_jets_p4 != nullptr){
                 for (size_t n = 0; n < other_jets_p4->size(); ++n){
                     LorentzVector1 other_jet = other_jets_p4->at(n);
                     unc->setJetPt(other_jet.pt());
                     unc->setJetEta(other_jet.eta());
                     const auto unc_result = unc->getUncertainty(scales.at(scale));
                     const double unc_var = unc_result ? *unc_result : 0;
                     const int sign = scales_variation.at(scale);
                     const auto sf = static_cast<typename LorentzVector1::Scalar>(1.0 + (sign * unc_var));
                     const auto shiftedMomentum = other_jet * sf;
                     shifted_met_px += other_jet.px() - shiftedMomentum.px();
                     shifted_met_py += other_jet.py() - shiftedMomentum.py();
                 }
            }

             shifted_met_px += met->px();
             shifted_met_py += met->py();
             double E = std::hypot(shifted_met_px,shifted_met_py);
             met->SetPxPyPzE(shifted_met_px,shifted_met_py,0,E);
         }

        return corrected_jets;
    }
private:
    std::map<analysis::UncertaintySource, std::shared_ptr<JetCorrectionUncertainty>> uncertainty_map;

};

} // namespace jec
