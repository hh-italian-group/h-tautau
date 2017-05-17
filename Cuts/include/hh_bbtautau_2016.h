/*! Double-Higgs to bbtautau selection.
Defined only cuts that are different from the H to tautau baseline selection.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace hh_bbtautau_2016 {

constexpr double DeltaR_betweenSignalObjects = 0.1; // > Decreased to not loose efficiency for high mX.

namespace MuTau {
    namespace muonID {
        constexpr double pt = 23; // > Increased to be away from the trigger threshold.
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
        constexpr double pt = 27; // > Increased to be away from the trigger threshold.

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
    constexpr double peak_tautau = 116;
    constexpr double resolution_tautau = 35.;
    constexpr double peak_bb = 111;
    constexpr double resolution_bb = 45;

    inline bool IsInsideEllipse (double mass_tautau, double mass_bb){
        const double ellipse_cut = std::pow(mass_tautau-peak_tautau, 2)/std::pow(resolution_tautau, 2) +
                std::pow(mass_bb-peak_bb, 2)/std::pow(resolution_bb, 2);
        return ellipse_cut<1;
    }
}

} // namespace hh_bbtautau_2016
} // namespace cuts
