/*! Definiton of HHGenEvent.
 This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/GenParticle.h"

namespace analysis {
    
    enum class TauGenDecayMode { Electron, Muon, Hadrons };

    ENUM_NAMES(TauGenDecayMode) = {
        { TauGenDecayMode::Electron, "Electron" },
        { TauGenDecayMode::Muon, "Muon" },
        { TauGenDecayMode::Hadrons, "Hadrons" }
    };
    
    struct HHGenEvent {
        const GenParticle *h_tautau, *h_bb;
        std::array<const GenParticle*, 2> tau, b;
        std::array<LorentzVectorM, 2> vis_tau;
        std::array<LorentzVectorXYZ, 2> vis_charged_tau;
        std::array<LorentzVectorXYZ, 2> b_jets;
        std::array<LorentzVectorXYZ, 5> b_jets_others;
        std::array<TauGenDecayMode, 2> tau_decay;
        LorentzVectorXYZ h_tautau_vis, h_bb_vis, h_bb_vis_all, h_bb_others_vis; //sum vis tau and bjets
    };
    
    template <class Vector1, class Vector2>
    inline bool hasDeltaRMatch(const Vector1& v1, const Vector2& v2, double deltaR)
    {
        return ROOT::Math::VectorUtil::DeltaR(v1, v2) < deltaR;

    }
}
