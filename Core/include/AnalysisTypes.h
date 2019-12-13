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

enum class LegType { e = 0, mu = 1, tau = 2, jet = 3 };
ENUM_NAMES(LegType) = {
    { LegType::e, "e" }, { LegType::mu, "mu" }, { LegType::tau, "tau" }, { LegType::jet, "jet" }
};

using ChannelLegTypes = std::pair<LegType, LegType>;
inline const ChannelLegTypes GetChannelLegTypes(Channel channel)
{
    static const std::map<Channel, ChannelLegTypes> leg_types {
        { Channel::ETau, { LegType::e, LegType::tau } }, { Channel::MuTau, { LegType::mu, LegType::tau } },
        { Channel::TauTau, { LegType::tau, LegType::tau } }, { Channel::MuMu, { LegType::mu, LegType::mu } },
    };
    auto iter = leg_types.find(channel);
    if(iter == leg_types.end())
        throw exception("Leg types are not defined for channel '%1%'") % channel;
    return iter->second;
}

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, TopPtUp = 5, TopPtDown = 6 };
ENUM_NAMES(EventEnergyScale) = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::TopPtUp, "TopPtUp" }, { EventEnergyScale::TopPtDown, "TopPtDown" },
};

enum class UncertaintyScale { Central = 0, Up = 1, Down = -1 };
ENUM_NAMES(UncertaintyScale) = {
    { UncertaintyScale::Central, "Central" }, { UncertaintyScale::Up, "Up" }, { UncertaintyScale::Down, "Down" }
};

enum class UncertaintySource { None = 0, TauES = 1, Total = 2, TopPt = 3, AbsoluteStat = 4, AbsoluteScale = 5,
     AbsoluteMPFBias = 6, AbsoluteFlavMap = 7, Fragmentation = 8, SinglePionECAL = 9, SinglePionHCAL = 10,
     FlavorQCD = 11, FlavorZJet = 12, FlavorPhotonJet = 13, FlavorPureGluon = 14, FlavorPureQuark = 15, FlavorPureCharm = 16,
     FlavorPureBottom = 17, TimePtEta = 18, RelativeJEREC1 = 19, RelativeJEREC2 = 20, RelativeJERHF = 21,
     RelativePtBB = 22, RelativePtEC1 = 23, RelativePtEC2 = 24, RelativePtHF = 25, RelativeBal = 26,
     RelativeFSR = 27, PileUpDataMC = 28, PileUpPtRef = 29, PileUpPtBB = 30, PileUpPtEC1 = 31, PileUpPtEC2 = 32,
     PileUpPtHF = 33, SubTotalPileUp = 34, SubTotalRelative = 35, SubTotalPt = 36, SubTotalScale = 37,
     SubTotalAbsolute = 38, SubTotalMC = 39, TotalNoFlavor = 40, TotalNoTime = 41, TotalNoFlavorNoTime = 42,
     EleFakingTauES = 43 };
ENUM_NAMES(UncertaintySource) = {
    { UncertaintySource::None, "None" }, { UncertaintySource::TauES, "TauES" },
    { UncertaintySource::Total, "Total" }, { UncertaintySource::TopPt, "TopPt" },
    { UncertaintySource::AbsoluteStat, "AbsoluteStat" },
    { UncertaintySource::AbsoluteScale, "AbsoluteScale" },
    { UncertaintySource::AbsoluteMPFBias, "AbsoluteMPFBias" },
    { UncertaintySource::AbsoluteFlavMap, "AbsoluteFlavMap" },
    { UncertaintySource::Fragmentation, "Fragmentation" }, { UncertaintySource::SinglePionECAL, "SinglePionECAL" },
    { UncertaintySource::SinglePionHCAL, "SinglePionHCAL" },
    { UncertaintySource::FlavorQCD, "FlavorQCD" },
    { UncertaintySource::FlavorZJet, "FlavorZJet" },
    { UncertaintySource::FlavorPhotonJet, "FlavorPhotonJet" },
    { UncertaintySource::FlavorPureGluon, "FlavorPureGluon" },
    { UncertaintySource::FlavorPureQuark, "FlavorPureQuark" },
    { UncertaintySource::FlavorPureCharm, "FlavorPureCharm" },
    { UncertaintySource::FlavorPureBottom, "FlavorPureBottom" },
    { UncertaintySource::TimePtEta, "TimePtEta" },
    { UncertaintySource::RelativeJEREC1, "RelativeJEREC1" },
    { UncertaintySource::RelativeJEREC2, "RelativeJEREC2" },
    { UncertaintySource::RelativeJERHF, "RelativeJERHF" },
    { UncertaintySource::RelativePtBB, "RelativePtBB" },
    { UncertaintySource::RelativePtEC1, "RelativePtEC1" },
    { UncertaintySource::RelativePtEC2, "RelativePtEC2" },
    { UncertaintySource::RelativePtHF, "RelativePtHF" },
    { UncertaintySource::RelativeBal, "RelativeBal" },
    { UncertaintySource::RelativeFSR, "RelativeFSR" },
    { UncertaintySource::PileUpDataMC, "PileUpDataMC" },
    { UncertaintySource::PileUpPtRef, "PileUpPtRef" },
    { UncertaintySource::PileUpPtBB, "PileUpPtBB" },
    { UncertaintySource::PileUpPtEC1, "PileUpPtEC1" },
    { UncertaintySource::PileUpPtEC2, "PileUpPtEC2" },
    { UncertaintySource::PileUpPtHF, "PileUpPtHF" },
    { UncertaintySource::SubTotalPileUp, "SubTotalPileUp" },
    { UncertaintySource::SubTotalRelative, "SubTotalRelative" },
    { UncertaintySource::SubTotalPt, "SubTotalPt" },
    { UncertaintySource::SubTotalScale, "SubTotalScale" },
    { UncertaintySource::SubTotalAbsolute, "SubTotalAbsolute" },
    { UncertaintySource::SubTotalMC, "SubTotalMC" },
    { UncertaintySource::TotalNoFlavor, "TotalNoFlavor"},
    { UncertaintySource::TotalNoTime, "TotalNoTime"},
    { UncertaintySource::TotalNoFlavorNoTime, "TotalNoFlavorNoTime"},
    { UncertaintySource::EleFakingTauES, "EleFakingTauES"}
};

using EventEnergyScaleSet = EnumNameMap<EventEnergyScale>::EnumEntrySet;

enum class DiscriminatorWP { VVVLoose = 0, VVLoose = 1, VLoose = 2, Loose = 3, Medium = 4, Tight = 5,
                             VTight = 6, VVTight = 7, VVVTight = 8 };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VVVLoose, "VVVLoose" }, { DiscriminatorWP::VVLoose, "VVLoose" },
    { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" }, { DiscriminatorWP::Medium, "Medium" },
    { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" }, { DiscriminatorWP::VVTight, "VVTight" },
    { DiscriminatorWP::VVVTight, "VVVTight" }
};
const EnumNameMap<DiscriminatorWP> __DiscriminatorWP_short_names("ShortWPNames", {
    { DiscriminatorWP::VVVLoose, "VVVL" }, { DiscriminatorWP::VVLoose, "VVL" }, { DiscriminatorWP::VLoose, "VL" },
    { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" }, { DiscriminatorWP::Tight, "T" },
    { DiscriminatorWP::VTight, "VT" }, { DiscriminatorWP::VVTight, "VVT" }, { DiscriminatorWP::VVVTight, "VVVT" }
});

enum class MetType { PF, MVA, PUPPI };
ENUM_NAMES(MetType) = {
    { MetType::PF, "PF" }, { MetType::MVA, "MVA" }, { MetType::PUPPI, "PUPPI" }
};

enum class Period { Run2015, Run2016, Run2017, Run2018 };
ENUM_NAMES(Period) = {
    { Period::Run2015, "Run2015" },
    { Period::Run2016, "Run2016" },
    { Period::Run2017, "Run2017" },
    { Period::Run2018, "Run2018" },
};

enum class GenLeptonMatch { Electron = 1, Muon = 2, TauElectron = 3,  TauMuon = 4, Tau = 5, NoMatch = 6 };
ENUM_NAMES(GenLeptonMatch) = {
    { GenLeptonMatch::Electron, "gen_electron" },
    { GenLeptonMatch::Muon, "gen_muon" },
    { GenLeptonMatch::TauElectron, "gen_electron_from_tau" },
    { GenLeptonMatch::TauMuon, "gen_muon_from_tau" },
    { GenLeptonMatch::Tau, "gen_tau" },
    { GenLeptonMatch::NoMatch, "no_gen_match" }
};

enum class GenEventType { Other = 0, TTbar_Hadronic = 1, TTbar_SemiLeptonic = 2, TTbar_Leptonic = 3 };
ENUM_NAMES(GenEventType) = {
    { GenEventType::Other, "other" },
    { GenEventType::TTbar_Hadronic, "TTbar_Hadronic" },
    { GenEventType::TTbar_SemiLeptonic, "TTbar_SemiLeptonic" },
    { GenEventType::TTbar_Leptonic, "TTbar_Leptonic" },
};

} // namespace analysis
