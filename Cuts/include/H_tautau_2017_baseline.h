/*! Higgs to tautau baseline selection for 2016 analyses.
If not specified otherwise, all definitions are taken from the Higgs2Tau Working TWiki for 2016:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace H_tautau_2017 {

constexpr double DeltaR_betweenSignalObjects = 0.5; // >
constexpr double DeltaR_triggerMatch = 0.5; // <

namespace MuTau {
    namespace muonID {
        constexpr double pt = 21; // >

    }

    namespace tauID {

    }

}

namespace ETau {
    namespace electronID {

    }

    namespace tauID {

    }


}

namespace TauTau {
    namespace tauID {

    }
}

namespace MuMu {
    namespace muonID {

    }

}



} // namespace H_tautau_2017
} // namespace cuts
