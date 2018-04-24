#include <string>
#include <vector>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "h-tautau/Analysis/include/Candidate.h"
#include "h-tautau/Analysis/include/TupleObjects.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

#include "DataFormats/METReco/interface/CorrMETData.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"

#include "JetMETCorrections/Type1MET/interface/AddCorrectionsToGenericMET.h"
#include "RecoMET/METAlgorithms/interface/METSignificance.h"

namespace jec {

    using JetCandidate = analysis::Candidate<ntuple::TupleJet>;
    using JetCollection = std::vector<JetCandidate>;
    using analysis::UncertaintySource;
    using analysis::UncertaintyScale;

    class JECUncertaintiesWrapper
    {
    public:


        static const std::set<UncertaintySource>& JetUncertainties()
        {
            static const std::set<UncertaintySource> jetUncertainties = {
                UncertaintySource::JetAbsolute, UncertaintySource::JetHighPtExtra,
                UncertaintySource::JetSinglePionECAL, UncertaintySource::JetSinglePionHCAL,
                UncertaintySource::JetFlavorQCD, UncertaintySource::JetTime,
                UncertaintySource::JetRelativeJEREC1, UncertaintySource::JetRelativeJEREC2,
                UncertaintySource::JetRelativeJERHF, UncertaintySource::JetRelativePtBB,
                UncertaintySource::JetRelativePtEC1, UncertaintySource::JetRelativePtEC2,
                UncertaintySource::JetRelativePtHF, UncertaintySource::JetRelativeFSR,
                UncertaintySource::JetRelativeStatEC2, UncertaintySource::JetRelativeStatHF,
                UncertaintySource::JetPileUpDataMC, UncertaintySource::JetPileUpPtBB,
                UncertaintySource::JetPileUpPtEC, UncertaintySource::JetPileUpPtHF,
                UncertaintySource::JetPileUpBias, UncertaintySource::JetSubTotalPileUp,
                UncertaintySource::JetSubTotalRelative, UncertaintySource::JetSubTotalPt,
                UncertaintySource::JetSubTotalMC, UncertaintySource::JetTotalNoFlavor,
                UncertaintySource::JetFlavorZJet, UncertaintySource::JetFlavorPhotonJet,
                UncertaintySource::JetFlavorPureGluon, UncertaintySource::JetFlavorPureQuark,
                UncertaintySource::JetFlavorPureCharm, UncertaintySource::JetFlavorPureBottom
            };
            return jetUncertainties;
        }

        static const std::set<UncertaintySource>& JetUncertainties_withTotal()
        {
            auto createUncSet = []() {
                std::set<UncertaintySource> jetUncertainties = JetUncertainties();
                jetUncertainties.insert(UncertaintySource::JetTotal);
                return jetUncertainties;
             };

            static const std::set<UncertaintySource> jetUncertaintiesTotal = createUncSet();
            return jetUncertaintiesTotal;
        }

        JECUncertaintiesWrapper(const std::string& uncertainties_source) {

           for (const auto jet_unc : JetUncertainties_withTotal()) {
               std::string full_name = analysis::ToString(jet_unc);
               const std::string name = full_name.substr(3);
               JetCorrectorParameters p(uncertainties_source, name);
               auto unc = std::make_shared<JetCorrectionUncertainty>(p);

               uncertainty_map[jet_unc] = unc;
           } // for unc

        }

        template<typename JetCollection, typename LorentzVector1, typename LorentzVector2>
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
                unc->setJetPt(jet.GetMomentum().pt());
                unc->setJetEta(jet.GetMomentum().eta());
                const double unc_var = unc->getUncertainty(scales.at(scale));
                const int sign = scales_variation.at(scale);
                const double sf = 1.0 + (sign * unc_var);
                const auto shiftedMomentum = jet.GetMomentum() * sf;
                JetCandidate corr_jet(jet);
                corr_jet.SetMomentum(shiftedMomentum);
                corrected_jets.push_back(corr_jet);
                shifted_met_px += jet.GetMomentum().px() - corr_jet.GetMomentum().px();
                shifted_met_py += jet.GetMomentum().py() - corr_jet.GetMomentum().py();
            }

            for (const auto other_jet : other_jets_p4){
                shifted_met_px -= other_jet->px();
                shifted_met_py -= other_jet->py();
            }

            if(met){
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

}
