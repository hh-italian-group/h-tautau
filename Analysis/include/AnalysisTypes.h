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

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4 };
ENUM_NAMES(EventEnergyScale) = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
};

//using EventEnergyScaleSet = EnumNameMap<EventEnergyScale>::EnumEntrySet;
//static const auto& AllEventEnergyScales = __EventEnergyScale_names<>::names.GetEnumEntries();

enum class DiscriminatorWP { VLoose, Loose, Medium, Tight, VTight };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" }, { DiscriminatorWP::Medium, "Medium" },
    { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" }
};

enum class MetType { PF, MVA, PUPPI };
ENUM_NAMES(MetType) = {
    { MetType::PF, "PF" }, { MetType::MVA, "MVA" }, { MetType::PUPPI, "PUPPI" }
};

enum class Period { Run2015, Run2016 };
ENUM_NAMES(Period) = {
    { Period::Run2015, "Run2015" },
    { Period::Run2016, "Run2016" }
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

} // namespace analysis
