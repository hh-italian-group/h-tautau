/*! Double-Higgs to bbtautau selection for the legacy Run 2 analysis.
Defined only cuts that are different from the H to tautau baseline selection.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "eleID_Run2.h"
#include "muonID_Run2.h"
#include "tauID_Run2.h"
#include "H_tautau_Run2.h"

namespace cuts {
namespace hh_bbtautau_Run2 {

constexpr double DeltaR_Lep_Lep = ::cuts::H_tautau_Run2::DeltaR_Lep_Lep; // >
constexpr double DeltaR_Lep_Jet = ::cuts::H_tautau_Run2::DeltaR_Lep_Jet; // >
constexpr double DeltaR_triggerMatch = ::cuts::H_tautau_Run2::DeltaR_triggerMatch; // <

namespace TauTau {
    namespace tauID {
        constexpr double pt = ::cuts::tauID_Run2::pt; // >
        constexpr double eta = ::cuts::tauID_Run2::eta_diTauTrigger; // <
        constexpr double dz = ::cuts::tauID_Run2::dz; // <
    }
}

namespace MuTau {
    namespace muonID {
        constexpr double pt = 18; // > [GeV] -- lowest pt_thr - 1 GeV
        constexpr double eta = ::cuts::H_tautau_Run2::MuTau::muonID::eta; // <
        constexpr double pfRelIso04 = ::cuts::H_tautau_Run2::MuTau::muonID::pfRelIso04; // <
        constexpr double dz = ::cuts::H_tautau_Run2::MuTau::muonID::dz; // <
        constexpr double dxy = ::cuts::H_tautau_Run2::MuTau::muonID::dxy; // <

        // pass Tight muon ID
    }

    namespace tauID {
        constexpr double pt = ::cuts::tauID_Run2::pt; // >
        constexpr double eta = ::cuts::tauID_Run2::eta; // <
        constexpr double dz = ::cuts::hh_bbtautau_Run2::TauTau::tauID::dz; // <
    }
}

namespace ETau {
    // ElectronID recommendation
    // https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_for_regul
    namespace electronID {
        constexpr double pt = 23; // > [GeV] -- lowest pt_thr - 1 GeV
        constexpr double eta = ::cuts::H_tautau_Run2::ETau::electronID::eta; // <
        constexpr double dz = ::cuts::H_tautau_Run2::ETau::electronID::dz; // <
        constexpr double dxy = ::cuts::H_tautau_Run2::ETau::electronID::dxy; // <
    }

    namespace tauID {
        constexpr double pt = ::cuts::hh_bbtautau_Run2::MuTau::tauID::pt; // >
        constexpr double eta = ::cuts::hh_bbtautau_Run2::MuTau::tauID::eta; // <
        constexpr double dz = ::cuts::hh_bbtautau_Run2::MuTau::tauID::dz; // <
    }
}

namespace MuMu {
    namespace muonID {
        constexpr double pt_sel = std::min(::cuts::hh_bbtautau_Run2::MuTau::muonID::pt,
                                           ::cuts::hh_bbtautau_Run2::MuTau::tauID::pt); // >
        constexpr double pt = 20; // > [GeV]
        constexpr double eta = std::max(::cuts::hh_bbtautau_Run2::MuTau::muonID::eta,
                                        ::cuts::hh_bbtautau_Run2::MuTau::tauID::eta); // <
        constexpr double dz = ::cuts::hh_bbtautau_Run2::MuTau::muonID::dz; // <
        constexpr double dxy = ::cuts::hh_bbtautau_Run2::MuTau::muonID::dxy; // <

        // pass Tight muon ID

        constexpr double pfRelIso04 = ::cuts::hh_bbtautau_Run2::MuTau::muonID::pfRelIso04; // <
        constexpr double m_vis_low = 60; // >
        constexpr double m_vis_high = 120; // <
    }
}

namespace electronVeto {
    constexpr double pt = ::cuts::H_tautau_Run2::electronVeto::pt; // >
    constexpr double eta = ::cuts::H_tautau_Run2::electronVeto::eta; // <
    constexpr double dz = ::cuts::H_tautau_Run2::electronVeto::dz; // <
    constexpr double dxy = ::cuts::H_tautau_Run2::electronVeto::dxy; // <

    // pass mvaEleID_iso_Medium OR (mvaEleID_noIso_Medium AND pfRelIso04 < 0.3)

    // missingHits and passConversionVeto should not be applied, because they are used as inputs to the MVA id and
    // direct cut on them is not recommented by POG.
}

namespace muonVeto {
    constexpr double pt = ::cuts::H_tautau_Run2::muonVeto::pt; // >
    constexpr double eta = ::cuts::H_tautau_Run2::muonVeto::eta; // <
    constexpr double dz = ::cuts::H_tautau_Run2::muonVeto::dz; // <
    constexpr double dxy = ::cuts::H_tautau_Run2::muonVeto::dxy; // <
    constexpr double pfRelIso04 = ::cuts::H_tautau_Run2::muonVeto::pfRelIso04; // <
    constexpr double tightIso = 0.1; // < (for pre-selection at the tuple production stage)

    // pass Medium OR Tight muon ID
}

namespace hh_tag {
    constexpr double boosted_m_tautau_min = 128;
    constexpr double boosted_m_tautau_max = 60;
    constexpr double boosted_m_bb_min = 159;
    constexpr double boosted_m_bb_max = 94;

    constexpr analysis::EllipseParameters m_hh_window{129.0, 53.0, 169.0, 145.0};
    constexpr analysis::EllipseParameters new_m_hh_window{87.9563, 41.8451, 109.639, 43.0346};

    inline constexpr bool IsInsideBoostedMassWindow(double mass_tautau, double mass_bb)
    {
        return mass_tautau > boosted_m_tautau_min && mass_tautau < boosted_m_tautau_max
            && mass_bb > boosted_m_bb_min && mass_bb < boosted_m_bb_max;
    }
}

namespace VBF {
    constexpr double mass_jj = 500;
    constexpr double deltaeta_jj = 3;
    constexpr double mass_jj_tight = 800;
}

namespace fatJetID {
    constexpr double mass = 30;
    constexpr double deltaR_subjet = 0.4;
}

namespace jetID {
    constexpr double pt = 20; // > [GeV]
    constexpr double pt_presel = 15; // > [GeV]
    constexpr double eta = 5; // <
    constexpr double vbf_pt = 30; // > [GeV]
    constexpr double vbf_eta = 4.7; // <

    //MET Recomendations. See presentation https://indico.cern.ch/event/762187/contributions/3218088/attachments/1753676/2842443/HTT_15Nov18_PileupJetID_ADow.pdf
    //and https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
    constexpr double max_pt_veto = 50; // < [GeV]
    constexpr double eta_low_veto = 2.65; // >
    constexpr double eta_high_veto = 3.139; // <
}

} // namespace hh_bbtautau_Run2
} // namespace cuts
