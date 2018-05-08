#include <string>
#include <vector>
#include "h-tautau/McCorrections/include/JetCorrectorParameters.h"
#include "h-tautau/McCorrections/include/JetCorrectionUncertainty.h"
#include "h-tautau/Analysis/include/Candidate.h"
#include "h-tautau/Analysis/include/TupleObjects.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

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
                UncertaintySource::AbsoluteStat,
                UncertaintySource::AbsoluteScale,
                UncertaintySource::AbsoluteMPFBias,
                UncertaintySource::AbsoluteFlavMap,
                UncertaintySource::Fragmentation,
                UncertaintySource::SinglePionECAL,
                UncertaintySource::SinglePionHCAL,
                UncertaintySource::FlavorQCD,
                UncertaintySource::FlavorZJet,
                UncertaintySource::FlavorPhotonJet,
                UncertaintySource::FlavorPureGluon,
                UncertaintySource::FlavorPureQuark,
                UncertaintySource::FlavorPureCharm,
                UncertaintySource::FlavorPureBottom,
                UncertaintySource::TimePtEta,
                UncertaintySource::RelativeJEREC1,
                UncertaintySource::RelativeJEREC2,
                UncertaintySource::RelativeJERHF,
                UncertaintySource::RelativePtBB,
                UncertaintySource::RelativePtEC1,
                UncertaintySource::RelativePtEC2,
                UncertaintySource::RelativePtHF,
                UncertaintySource::RelativeBal,
                UncertaintySource::RelativeFSR,
                UncertaintySource::PileUpDataMC,
                UncertaintySource::PileUpPtRef,
                UncertaintySource::PileUpPtBB,
                UncertaintySource::PileUpPtEC1,
                UncertaintySource::PileUpPtEC2,
                UncertaintySource::PileUpPtHF,
                UncertaintySource::SubTotalPileUp,
                UncertaintySource::SubTotalRelative,
                UncertaintySource::SubTotalPt,
                UncertaintySource::SubTotalScale,
                UncertaintySource::SubTotalAbsolute,
                UncertaintySource::SubTotalMC,
                UncertaintySource::TotalNoFlavor,
                UncertaintySource::TotalNoTime,
                UncertaintySource::TotalNoFlavorNoTime
            };
            return jetUncertainties;
        }

        static const std::set<UncertaintySource>& JetUncertainties_withTotal()
        {
            auto createUncSet = []() {
                std::set<UncertaintySource> jetUncertainties = JetUncertainties();
                jetUncertainties.insert(UncertaintySource::Total);
                return jetUncertainties;
             };

            static const std::set<UncertaintySource> jetUncertaintiesTotal = createUncSet();
            return jetUncertaintiesTotal;
        }

        JECUncertaintiesWrapper(const std::string& uncertainties_source) {

            for (const auto jet_unc : JetUncertainties_withTotal()) {
                std::string full_name = analysis::ToString(jet_unc);
//                const std::string name = full_name.substr(3);
                JetCorrectorParameters p(uncertainties_source, full_name);
                auto unc = std::make_shared<JetCorrectionUncertainty>(p);
           
                uncertainty_map[jet_unc] = unc;
            } // for unc

        }

        template<typename JetCollection, typename LorentzVector1 = analysis::LorentzVector, typename LorentzVector2 = analysis::LorentzVector>
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
            
             const double unc_var = unc->getUncertainty(scales.at(scale));
             const int sign = scales_variation.at(scale);
            
             for (const auto& jet : jet_candidates){
                 unc->setJetPt(static_cast<float>(jet.GetMomentum().pt()));
                 unc->setJetEta(static_cast<float>(jet.GetMomentum().eta()));
                 const double sf = 1.0 + (sign * unc_var);
                 const auto shiftedMomentum = jet.GetMomentum() * sf;
                 JetCandidate corr_jet(jet);
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
                         const double unc_var = unc->getUncertainty(scales.at(scale));
                         const int sign = scales_variation.at(scale);
                         const double sf = 1.0 + (sign * unc_var);
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

}
