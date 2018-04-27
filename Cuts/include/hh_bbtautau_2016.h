/*! Double-Higgs to bbtautau selection.
Defined only cuts that are different from the H to tautau baseline selection.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace cuts {
namespace hh_bbtautau_2016 {

constexpr double DeltaR_betweenSignalObjects = 0.1; // > Decreased to not loose efficiency for high mX.

namespace MuTau {
    namespace muonID {
        //constexpr double pt = 23; // > Increased to be away from the trigger threshold.
        constexpr double pt = 19; // > Lowest possible trigger threshold plus safetyPt.
        constexpr bool isTightMuon = true; // = Should be ok, since the HIP problem is solved.

        // pfRelIso should be applied at the tuple production level.
    }

    namespace tauID {
        // againstElectronVLooseMVA6 and againstMuonTight3 should be applied at the tuple production level.
    }

    // ZmumuVeto should not be applied.
}

namespace ETau {
    namespace electronID {
        //constexpr double pt = 27; // > Increased to be away from the trigger threshold.
        constexpr double pt = 25; // > Lowest possible trigger threshold plus safetyPt.

        // missingHits and passConversionVeto should not be applied, because they are used as inputs to the MVA id and
        // direct cut on them is not recommented by POG.

        // pfRelIso should be applied at the tuple production level.
    }

    namespace tauID {
        // againstElectronTightMVA6 and againstMuonLoose3 should be applied at the tuple production level.
    }

    // ZeeVeto should not be applied.
}

namespace MuMu {
    namespace muonID {
        constexpr double pt = 10; // >
        constexpr bool isTightMuon = true; // = Same as for muTau channel
        constexpr double eta_leading = 2.1; // < Same as for muTau channel

        // pfRelIso should be applied at the tuple production level.
    }

    constexpr double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects; // >
}

namespace electronVeto {
    // missingHits and passConversionVeto should not be applied, because they are used as inputs to the MVA id and
    // direct cut on them is not recommented by POG.
}

namespace muonVeto {
    constexpr bool isLooseMuon = true; // = Should be ok, since the HIP problem is solved.
}

namespace hh_tag {

    constexpr double boosted_m_tautau_min = 80;
    constexpr double boosted_m_tautau_max = 152;
    constexpr double boosted_m_bb_min = 90;
    constexpr double boosted_m_bb_max = 160;

    inline const analysis::EllipseParameters& m_hh_window()
    {
        static const analysis::EllipseParameters m_hh_window_param{116.0, 35.0, 111.0, 45.0};
        return m_hh_window_param;
    }

    inline const analysis::EllipseParameters& new_m_hh_window()
    {
        static const analysis::EllipseParameters m_hh_window_param{87.9563, 41.8451, 109.639, 43.0346};
        return m_hh_window_param;
    }

    inline bool IsInsideBoostedMassWindow(double mass_tautau, double mass_bb)
    {
        return mass_tautau > boosted_m_tautau_min && mass_tautau < boosted_m_tautau_max
            && mass_bb > boosted_m_bb_min && mass_bb < boosted_m_bb_max;
    }
}

namespace fatJetID {
    constexpr double mass = 30;
    constexpr double deltaR_subjet = 0.4;
}

} // namespace hh_bbtautau_2016
} // namespace cuts
