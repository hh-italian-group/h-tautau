/*! Common simple types for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2, MuMu = 3 };
ENUM_NAMES(Channel) = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }, { Channel::MuMu, "muMu" }
};
const EnumNameMap<Channel> __Channel_names_latex("ChannelLatex", {
    { Channel::ETau, "e#tau_{h}" }, { Channel::MuTau, "#mu#tau_{h}" }, { Channel::TauTau, "#tau_{h}#tau_{h}" },
    { Channel::MuMu, "#mu#mu" }
});

enum class LegType { e = 0, mu = 1, tau = 2 };
ENUM_NAMES(LegType) = {
    { LegType::e, "e" }, { LegType::mu, "mu" }, { LegType::tau, "tau" }
};

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, TopPtUp = 5, TopPtDown = 6 };
ENUM_NAMES(EventEnergyScale) = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::TopPtUp, "TopPtUp" }, { EventEnergyScale::TopPtDown, "TopPtDown" },
};

enum class UncertaintyScale { Central = 0, Up = 1, Down = 2 };
ENUM_NAMES(UncertaintyScale) = {
    { UncertaintyScale::Central, "Central" }, { UncertaintyScale::Up, "Up" }, { UncertaintyScale::Down, "Down" }
};

enum class UncertaintySource { None = 0, TauES = 1, JetTotal = 2, TopPt = 3, JetAbsolute = 4,
    JetHighPtExtra = 5, JetSinglePionECAL = 6, JetSinglePionHCAL = 7, JetFlavorQCD = 8, JetTime = 9,
   JetRelativeJEREC1 = 10, JetRelativeJEREC2 = 11, JetRelativeJERHF = 12,
   JetRelativePtBB = 13, JetRelativePtEC1 = 14, JetRelativePtEC2 = 15, JetRelativePtHF = 16,
   JetRelativeFSR = 17, JetRelativeStatEC2 = 18, JetRelativeStatHF = 19, JetPileUpDataMC = 20,
   JetPileUpPtBB = 21, JetPileUpPtEC = 22, JetPileUpPtHF = 23, JetPileUpBias = 24,
   JetSubTotalPileUp = 25, JetSubTotalRelative = 26, JetSubTotalPt = 27, JetSubTotalMC = 28,
   JetTotalNoFlavor = 29, JetFlavorZJet = 30, JetFlavorPhotonJet = 31, JetFlavorPureGluon = 32,
   JetFlavorPureQuark = 33, JetFlavorPureCharm = 34, JetFlavorPureBottom = 35 };
ENUM_NAMES(UncertaintySource) = {
    { UncertaintySource::None, "None" }, { UncertaintySource::TauES, "TauES" },
    { UncertaintySource::JetTotal, "JetTotal" }, { UncertaintySource::TopPt, "TopPt" },
    { UncertaintySource::JetAbsolute, "JetAbsolute" },
    { UncertaintySource::JetHighPtExtra, "JetHighPtExtra" },
    { UncertaintySource::JetSinglePionECAL, "JetSinglePionECAL" },
    { UncertaintySource::JetSinglePionHCAL, "JetSinglePionHCAL" },
    { UncertaintySource::JetFlavorQCD, "JetFlavorQCD" }, { UncertaintySource::JetTime, "JetTime" },
    { UncertaintySource::JetRelativeJEREC1, "JetRelativeJEREC1" },
    { UncertaintySource::JetRelativeJEREC2, "JetRelativeJEREC2" },
    { UncertaintySource::JetRelativeJERHF, "JetRelativeJERHF" },
    { UncertaintySource::JetRelativePtBB, "JetRelativePtBB" },
    { UncertaintySource::JetRelativePtEC1, "JetRelativePtEC1" },
    { UncertaintySource::JetRelativePtEC2, "JetRelativePtEC2" },
    { UncertaintySource::JetRelativePtHF, "JetRelativePtHF" },
    { UncertaintySource::JetRelativeFSR, "JetRelativeFSR" },
    { UncertaintySource::JetRelativeStatEC2, "JetRelativeStatEC2" },
    { UncertaintySource::JetRelativeStatHF, "JetRelativeStatHF" },
    { UncertaintySource::JetPileUpDataMC, "JetPileUpDataMC" },
    { UncertaintySource::JetPileUpPtBB, "JetPileUpPtBB" },
    { UncertaintySource::JetPileUpPtEC, "JetPileUpPtEC" },
    { UncertaintySource::JetPileUpPtHF, "JetPileUpPtHF" },
    { UncertaintySource::JetPileUpBias, "JetPileUpBias" },
    { UncertaintySource::JetSubTotalPileUp, "JetSubTotalPileUp" },
    { UncertaintySource::JetSubTotalRelative, "JetSubTotalRelative" },
    { UncertaintySource::JetSubTotalPt, "JetSubTotalPt" },
    { UncertaintySource::JetSubTotalMC, "JetSubTotalMC" },
    { UncertaintySource::JetTotalNoFlavor, "JetTotalNoFlavor" },
    { UncertaintySource::JetFlavorZJet, "JetFlavorZJet" },
    { UncertaintySource::JetFlavorPhotonJet, "JetFlavorPhotonJet" },
    { UncertaintySource::JetFlavorPureGluon, "JetFlavorPureGluon" },
    { UncertaintySource::JetFlavorPureQuark, "JetFlavorPureQuark" },
    { UncertaintySource::JetFlavorPureCharm, "JetFlavorPureCharm" },
    { UncertaintySource::JetFlavorPureBottom, "JetFlavorPureBottom" }

};

using EventEnergyScaleSet = EnumNameMap<EventEnergyScale>::EnumEntrySet;

enum class DiscriminatorWP { VLoose, Loose, Medium, Tight, VTight, VVTight };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" }, { DiscriminatorWP::Medium, "Medium" },
    { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" }, { DiscriminatorWP::VVTight, "VVTight" }
};

enum class MetType { PF, MVA, PUPPI };
ENUM_NAMES(MetType) = {
    { MetType::PF, "PF" }, { MetType::MVA, "MVA" }, { MetType::PUPPI, "PUPPI" }
};

enum class Period { Run2015, Run2016, Run2017 };
ENUM_NAMES(Period) = {
    { Period::Run2015, "Run2015" },
    { Period::Run2016, "Run2016" },
    { Period::Run2017, "Run2017" }
};

enum class GenMatch { Electron = 1, Muon = 2, TauElectron = 3,  TauMuon = 4, Tau = 5, NoMatch = 6 };
ENUM_NAMES(GenMatch) = {
    { GenMatch::Electron, "gen_electron" },
    { GenMatch::Muon, "gen_muon" },
    { GenMatch::TauElectron, "gen_electron_from_tau" },
    { GenMatch::TauMuon, "gen_muon_from_tau" },
    { GenMatch::Tau, "gen_tau" },
    { GenMatch::NoMatch, "no_gen_match" }
};

enum class GenEventType { Other = 0, TTbar_Hadronic = 1, TTbar_SemiLeptonic = 2, TTbar_Leptonic = 3 };
ENUM_NAMES(GenEventType) = {
    { GenEventType::Other, "other" },
    { GenEventType::TTbar_Hadronic, "TTbar_Hadronic" },
    { GenEventType::TTbar_SemiLeptonic, "TTbar_SemiLeptonic" },
    { GenEventType::TTbar_Leptonic, "TTbar_Leptonic" },
};

} // namespace analysis
