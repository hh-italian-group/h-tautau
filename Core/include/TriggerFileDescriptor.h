/*! Definition of the file descriptor for Trigger setup and descriptor.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace trigger_tools {

struct SetupDescriptor {
    std::string name;
    std::map<analysis::LegType, double> deltaPt_map;

};

using SetupDescriptorCollection = std::unordered_map<std::string, SetupDescriptor>;

struct TriggerFileDescriptor {
    std::string name;
    std::set<analysis::Channel> channels;
    std::vector<std::string> legs;

};

using TriggerFileDescriptorCollection = std::unordered_map<std::string, TriggerFileDescriptor>;

} // namespace analysis
