/*! Higgs to tautau baseline selection for 2016 analyses.
If not specified otherwise, all definitions are taken from the Higgs2Tau Working TWiki for 2016:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace H_tautau_2016 {

constexpr double DeltaR_betweenSignalObjects = 0.5; // >
constexpr double DeltaR_triggerMatch = 0.5; // <

namespace MuTau {
    namespace muonID {
        constexpr double pt = 20; // >
        constexpr double eta = 2.1; // <
        constexpr double dz = 0.2; // <
        constexpr double dxy = 0.045; // <

        // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Short_Term_Medium_Muon_Definitio
        constexpr bool isMediumMuon = true; // =

        // Post-sync ntuple cuts

        constexpr double pfRelIso04 = 0.15; // <
    }

    namespace tauID {
        constexpr double pt = 20; // >
        constexpr double eta = 2.3; // <
        constexpr double decayModeFinding = 0.5; // >
        constexpr double dz = 0.2; // <
        constexpr int absCharge = 1; // =

        // Post-sync ntuple cuts

        constexpr double againstElectronVLooseMVA6 = 0.5; // >
        constexpr double againstMuonTight3 = 0.5; // >

        // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        constexpr double byTightIsolationMVArun2v1DBoldDMwLT = 0.5; // >
    }

    namespace ZmumuVeto {
        constexpr double pt = 15; // >
        constexpr double eta = 2.4; // <
        constexpr double dz = 0.2; // <
        constexpr double dxy = 0.045; // <
        constexpr bool isGlobalMuon = true; // =
        constexpr bool isTrackerMuon = true; // =
        constexpr bool isPFMuon = true; // =
        constexpr double pfRelIso04 = 0.3; // <
        constexpr double deltaR = 0.15; // >
        constexpr bool haveOppositeCharge = true; // =
    }
}

namespace ETau {
    namespace electronID {
        constexpr double pt = 25; // >
        constexpr double eta = 2.1; // <
        constexpr double dz = 0.2; // <
        constexpr double dxy = 0.045; // <

        // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
        constexpr bool MVApogID80effWP = true; // =

        constexpr int missingHits = 1; // <
        constexpr bool passConversionVeto = true; // =

        // Post-sync ntuple cuts
        constexpr double pfRelIso04 = 0.1; // <
    }

    namespace tauID {
        constexpr double pt = 20; // >
        constexpr double eta = 2.3; // <
        constexpr double decayModeFinding = 0.5; // >
        constexpr double dz = 0.2; // <
        constexpr int absCharge = 1; // =

        // Post-sync ntuple cuts

        constexpr double againstElectronTightMVA6 = 0.5; // >
        constexpr double againstMuonLoose3 = 0.5; // >

        // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        constexpr double byTightIsolationMVArun2v1DBoldDMwLT = 0.5; // >
    }

    namespace ZeeVeto {
        constexpr double pt = 15; // >
        constexpr double eta = 2.5; // <
        constexpr double dz = 0.2; // <
        constexpr double dxy = 0.045; // <
        constexpr double pfRelIso04 = 0.3; // <
        constexpr double deltaR = 0.15; // >
        constexpr bool haveOppositeCharge = true; // =

        // Defined in https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        constexpr bool POGcutBasedVetoID = true; // =
    }
}

namespace TauTau {
    namespace tauID {
        constexpr double pt = 40; // >
        constexpr double eta = 2.1; // <
        constexpr double dz = 0.2; // <
        constexpr double decayModeFinding = 0.5; // >
        constexpr int absCharge = 1; // =

        // Post-sync ntuple cuts

        constexpr double againstElectronVLooseMVA6 = 0.5; // >
        constexpr double againstMuonLoose3 = 0.5; // >

        // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        constexpr double byVTightIsolationMVArun2v1DBoldDMwLT = 0.5; // >
    }
}

namespace electronVeto {
    constexpr double pt = 10; // >
    constexpr double eta = 2.5; // <
    constexpr double dz = 0.2; // <
    constexpr double dxy = 0.045; // <
    constexpr double pfRelIso04 = 0.3; // <

    // Defined in https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
    constexpr bool MVApogID90effWP = true; // =

    constexpr bool passConversionVeto = true; // =
    constexpr int missingHits = 1; // <=
}

namespace muonVeto {
    constexpr double pt = 10; // >
    constexpr double eta = 2.4; // <
    constexpr double dz = 0.2; // <
    constexpr double dxy = 0.045; // <
    constexpr bool isMediumMuon = true; // = HIP safe
    constexpr double pfRelIso04 = 0.3; // <
}

namespace jetID {
    constexpr double pt = 20; // >
    constexpr double eta = 4.7; // <
    constexpr bool pfLooseID = true; // =
    constexpr double deltaR_signalObjects = 0.5; // >
}

namespace vertex {
    // Take the primary vertex as the first one in the vertex collection, as these have been pre-sorted.
}

} // namespace H_tautau_2016
} // namespace cuts
