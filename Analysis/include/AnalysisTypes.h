/*! Common simple types for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, BtagEfficiencyUp = 5,
                              BtagEfficiencyDown = 6 , BtagFakeUp = 7, BtagFakeDown = 8};

enum class DiscriminatorWP { VLoose, Loose, Medium, Tight, VTight };
enum class MetType { PF, MVA, PUPPI };

ENUM_NAMES(Channel) = { { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" } };
const EnumNameMap<Channel> __Channel_names_latex("ChannelLatex", {
    { Channel::ETau, "e#tau_{h}" }, { Channel::MuTau, "#mu#tau_{h}" }, { Channel::TauTau, "#tau_{h}#tau_{h}" }
});

ENUM_NAMES(EventEnergyScale) = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::BtagEfficiencyUp, "BtagEfficiencyUp" },
    { EventEnergyScale::BtagEfficiencyDown, "BtagEfficiencyDown" },
    { EventEnergyScale::BtagFakeUp, "BtagFakeUp" },
    { EventEnergyScale::BtagFakeDown, "BtagFakeDown" }
};

using EventEnergyScaleSet = EnumNameMap<EventEnergyScale>::EnumEntrySet;
static const auto& AllEventEnergyScales = __EventEnergyScale_names.GetEnumEntries();

} // namespace analysis
