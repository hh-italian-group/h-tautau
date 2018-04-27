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

using ChannelLegTypes = std::pair<LegType, LegType>;
inline const ChannelLegTypes GetChannelLegTypes(Channel channel)
{
    static const std::unordered_map<Channel, ChannelLegTypes> leg_types {
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

enum class UncertaintyScale { Central = 0, Up = 1, Down = 2 };
ENUM_NAMES(UncertaintyScale) = {
    { UncertaintyScale::Central, "Central" }, { UncertaintyScale::Up, "Up" }, { UncertaintyScale::Down, "Down" }
};

using EventEnergyScaleSet = EnumNameMap<EventEnergyScale>::EnumEntrySet;

enum class DiscriminatorWP { VVLoose, VLoose, Loose, Medium, Tight, VTight, VVTight };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VVLoose, "VVLoose" }, { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" },
    { DiscriminatorWP::Medium, "Medium" }, { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" },
    { DiscriminatorWP::VVTight, "VVTight" }
};
const EnumNameMap<DiscriminatorWP> __DiscriminatorWP_short_names("ShortWPNames", {
    { DiscriminatorWP::VVLoose, "VVL" }, { DiscriminatorWP::VLoose, "VL" }, { DiscriminatorWP::Loose, "L" },
    { DiscriminatorWP::Medium, "M" }, { DiscriminatorWP::Tight, "T" }, { DiscriminatorWP::VTight, "VT" },
    { DiscriminatorWP::VVTight, "VVT" }
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
