/*! Common simple types for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/PhysicalValue.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, BtagEfficiencyUp = 5,
                              BtagEfficiencyDown = 6 , BtagFakeUp = 7, BtagFakeDown = 8};

namespace detail {
const std::map<Channel, std::string> ChannelNameMap = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }
};

const std::map<Channel, std::string> ChannelNameMapLatex = {
    { Channel::ETau, "e#tau_{h}" }, { Channel::MuTau, "#mu#tau_{h}" }, { Channel::TauTau, "#tau_{h}#tau_{h}" }
};

const std::map<EventEnergyScale, std::string> EventEnergyScaleNameMap = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::BtagEfficiencyUp, "BtagEfficiencyUp" },
    { EventEnergyScale::BtagEfficiencyDown, "BtagEfficiencyDown" },
    { EventEnergyScale::BtagFakeUp, "BtagFakeUp" },
    { EventEnergyScale::BtagFakeDown, "BtagFakeDown" }
};

} // namespace detail

typedef std::set<EventEnergyScale> EventEnergyScaleSet;

const std::set<EventEnergyScale> AllEventEnergyScales = tools::collect_map_keys(detail::EventEnergyScaleNameMap);

std::ostream& operator<< (std::ostream& s, const Channel& c)
{
    s << detail::ChannelNameMap.at(c);
    return s;
}

std::istream& operator>> (std::istream& s, Channel& c)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::ChannelNameMap) {
        if(map_entry.second == name) {
            c = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown channel name '") << name << "'.";
}

std::ostream& operator<< (std::ostream& s, const EventEnergyScale& es)
{
    s << detail::EventEnergyScaleNameMap.at(es);
    return s;
}

std::istream& operator>> (std::istream& s, EventEnergyScale& es)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::EventEnergyScaleNameMap) {
        if(map_entry.second == name) {
            es = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown event energy scale '") << name << "'.";
}


} // namespace analysis
