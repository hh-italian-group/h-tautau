/*! Definition of the file descriptor for Trigger setup and descriptor.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include "AnalysisTools/Core/include/map_vec.h"
#include "h-tautau/Core/include/AnalysisTypes.h"

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
    bool apply_data{true}, apply_mc{true};
    boost::optional<unsigned> min_run, max_run;
};

using TriggerFileDescriptorCollection = analysis::map_vec<std::string, TriggerFileDescriptor>;

} // namespace analysis
