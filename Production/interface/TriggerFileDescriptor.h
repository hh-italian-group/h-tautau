/*! Definition of the file descriptor for DY and Wjets sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "AnalysisTools/Core/include/ConfigReader.h"

namespace analysis {

namespace trigger{

struct SetupDescriptor {
    std::string name;
    std::map<analysis::LegType, double> deltaPt_map;

};

using SetupDescriptorCollection = std::unordered_map<std::string, SetupDescriptor>;

struct TriggerFileDescriptor {
    std::string name;
    std::set<Channel> channels;
    std::vector<std::string> legs;

};

using TriggerFileDescriptorCollection = std::unordered_map<std::string, TriggerFileDescriptor>;

} //namespace trigger

} // namespace analysis
