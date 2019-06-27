/*! Definiton of HHGenEvent.
 This file is part of https://github.com/hh-italian-group/h-tautau. */

namespace analysis {
    
    enum class GenDecayMode { Electron, Muon, Hadrons };
    
    struct HHGenEvent {
        const GenParticle *h_tautau, *h_bb;
        std::array<const GenParticle*, 2> tau, b;
        std::array<LorentzVectorM, 2> vis_tau;
        std::array<LorentzVectorXYZ, 2> vis_charged_tau;
        std::array<LorentzVectorXYZ, 2> b_jets;
        std::array<LorentzVectorXYZ, 5> b_jets_all;
        std::array<GenDecayMode, 2> tau_decay;
        LorentzVectorXYZ h_tautau_vis, h_bb_vis, h_bb_vis_all; //sum vis tau and bjets
    };
    
    template <class Vector1, class Vector2>
    inline const bool HasMatchWithMCObject(const Vector1 v1, const Vector2 v2, double deltaR)
    {
        if(ROOT::Math::VectorUtil::DeltaR(v1, v2) < deltaR)
            return true;
        else
            return false;
    }
}
