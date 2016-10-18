/*! Higgs in tautau recommended baseline selection cuts for 2015.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include <string>
#include <vector>
#include <map>

#include <TMath.h>
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/exception.h"

namespace cuts {
namespace Htautau_2015 {
// AN-2013/178 H->etau,mutau
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/analyzers/DiLeptonAnalyzer.py
const double DeltaR_betweenSignalObjects = 0.0001; // >

// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauMuAnalyzer.py#L272
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L665
const double DeltaR_triggerMatch = 0.5; // <

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#Baseline_e_tau_h_AN1
const double DeltaR_DileptonVeto = 0.15;

namespace MuTau {
    namespace muonID {
        const double pt = 23; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.1; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dB = 0.045; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        const bool isMediumMuon = true; // = twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
                                        // def : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Medium_Muon

        //After Synctuple

        const double pFRelIso = 0.15; // < twiki HiggsToTauTauWorking2015#Muons

        const double mt = 30; // < not used for sync, only for the final selection.

    }

    namespace tauID {
        const double pt = 20; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.3; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double decayModeFinding = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        //After Synctuple

        const double againstMuonTight3 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronVLooseMVA5 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronVLooseMVA6 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h 76x
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

   }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
    namespace ZmumuVeto {
        const double pt = 15; // >
        const double eta = 2.4; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // < same definition as db
        const bool isGlobalMuon = true; // = already applied at Tree selection level as skim - MuonBlock.cc
        const bool isTrackerMuon = true; // =
        const bool isPFMuon = true; // =
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =
    }

}

namespace ETau {
    // twiki HiggsToTauTauWorkingSummer2013
    namespace electronIDscaleFactor {
        const std::vector<double> eta = { 1.479, 2.1 };
        const std::vector<double> pt = { 24, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.8999, 0.7945 },
                                                                    { 0.9486, 0.8866 } };
    }

    // twiki HiggsToTauTauWorkingSummer2013
    namespace electronISOscaleFactor {
        const std::vector<double> eta = { 1.479, 2.1 };
        const std::vector<double> pt = { 24, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.9417, 0.9471 },
                                                                    { 0.9804, 0.9900 } };
    }

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#Baseline_e_tau_h_AN1
    namespace electronID{
        const double pt = 27; // >
        const double eta = 2.1; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // <
        const bool MVApogNonTrig = true; // tight Use Value map compute by POG people
        const int missingHits = 1; // <
        const bool passConversionVeto = true; // =

        const double pFRelIso = 0.1; // < twiki
        const double scEta_min[2] = { 0.8, 1.479 }; // tight


        const double mt = 30; // < not used for sync, only for the final selection.
    }

    namespace tauID {
        const double pt = 20; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.3; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double decayModeFinding = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        //After Synctuple

        const double againstMuonLoose3 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronTightMVA6 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

   }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
    namespace ZeeVeto {
        const double pt = 15; // >
        const double eta = 2.5; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // <
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =

        // WP95 definition for Veto twiki EgammaCutBasedIdentification
        const double barrel_eta_high = 1.479; // <=
        const double endcap_eta_low = 1.479; // >
        const double endcap_eta_high = 2.5; // <
        const size_t barrel_index = 0, endcap_index = 1;
        const double sigma_ieta_ieta[] = { 0.01, 0.03 }; // <
        const double delta_eta[] = { 0.007, 0.01 }; // <
        const double delta_phi[] = { 0.8, 0.7 }; // <
        const double HoverE[] = { 0.15, std::numeric_limits<double>::max() }; // <
        const double dZ_vtx[] = { 0.2, 0.2 }; // <
    }

}

namespace TauTau {
    namespace tauID {
        const double pt = 40; // > AN-2013/189 H->tautau full hadronic
        const double eta = 2.1; // < AN-2013/189 H->tautau full hadronic
        const double dz = 0.2; // < AN-2013/189 H->tautau full hadronic
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA3 = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading pt tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more
                                                     // for new tauID
        //After Synctuple

        const double againstMuonTight3 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronVLooseMVA6 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double byMediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1; // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                //https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/RecoTauTag/Configuration/python/HPSPFTaus_cff.py#L183

        namespace BackgroundEstimation {
            const double Isolation_upperLimit = 4; // < upper value of isolation for QCD data driven estimation
        }

    }
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorking2015#Third_lepton_vetoes
namespace electronVeto {
    const double pt = 10; // >
    const double eta_high = 2.5; // <
    const double dz = 0.2; // <
    const double d0 = 0.045; // <
    const double pFRelIso = 0.3; // <
    const double ref_pt = 20; // twiki HiggsToTauTauWorkingSummer2013#Electron_ID
    const double scEta_min[2] = {0.8, 1.479}; // loose HiggsToTauTauWorkingSummer2013#Electron_ID
    const double MVApogNonTrig[2][3] = {{0.925, 0.915, 0.965},{0.905,0.955, 0.975}};
                                              // loose HiggsToTauTauWorkingSummer2013#Electron_ID
    const int missingHits = 1; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
    const bool hasMatchedConversion = false; // =  HiggsToTauTauWorkingSummer2013#Electron_ID
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py
namespace muonVeto {
    const double pt = 10; // >
    const double eta = 2.4; // <
    const double dz = 0.2; // <
    const double dB = 0.045; // <

    //def of isMediumMuon: twiki SWGuideMuonId#Medium_Muon
    const bool isMediumMuon = true; // =

    const double pfRelIso = 0.3; // <
}

namespace jetID {
    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 30; // > , njets is filled with a 30 GeV cuts. All the other variables use jets with 20 pt cut
    const double eta = 4.7; // <
    const bool puLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >

    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const bool pfLooseID = true; // =

    const double pt_loose = 20; // >

}

namespace btag {
    // pfCombinedInclusiveSecondaryVertexV2BJetTags, twiki
    const double CSVL = 0.605; // > loose
    const double CSVM = 0.89; // > medium
    const double CSVT = 0.97; // > tight

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 20; // >
    const double eta = 2.4; // <
    const double CSV = 0.8; // > HTauTau Twiki
    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const bool puLooseID = true; // =
    const bool pfLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >
}

// AN-2013/188 H->tautau physics objects
namespace vertex {
    const double ndf = 4; // >
    const double z = 24.0; // < cm
    const double r = 2.0; // < cm
    const bool chooseHighestSumPt2 = true; // =
}

namespace tauCorrections {

    const double DecayModeWeight = 0.88; // = HiggsToTauTauWorkingSummer2013#TauTau_scale_factors
                                         // for 1-prong no pi 0 taus

    const double deltaR_matchGenParticle = 0.5; // gen Particle match
    const double deltaR = 0.3; // < Updated to be compatible with H->tautau code

    const double energyUncertainty = 0.03;
}

namespace jetToTauFakeRateWeight {

    inline double CalculateJetToTauFakeWeight(double tauPt)
    {
        const double tau_pt = tauPt < 200 ? tauPt : 200 ;
        return (1.15743)-(0.00736136*tau_pt)+(4.3699e-05*tau_pt*tau_pt)-(1.188e-07*tau_pt*tau_pt*tau_pt);
    }
}

namespace DrellYannCategorization {
    const double minimal_genParticle_pt = 8; // > GeV
    const double deltaR_matchGenParticle = 0.5; // < twiki HiggsToTauTauWorkingSummer2013#E_MU_TAU_channel
    // < 0.3 https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L757

    const double minimal_visible_momentum = 18; // > GeV
}

namespace DYEmbedded {
    namespace trigger {
        // Riccardo repository: https://github.com/rmanzoni/HTT/tree/
        //a03b227073b2d4d8a2abe95367c014694588bf98/CMGTools/H2TauTau/python/proto/samples/run2012/trigger*.py
        const std::set<std::string> hltPaths = { "HLT_Mu17_Mu8"};
    }

    const double invariantMassCut = 50.; // > GeV
           // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/EmbedWeighter.py#L132

    const double deltaR_matchGenParticle = DrellYannCategorization::deltaR_matchGenParticle;
}

} // Htautau_2015
} // cuts
