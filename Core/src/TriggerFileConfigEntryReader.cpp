/*! Definition of the file configuration entry reader for Trigger configuration file.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "h-tautau/Core/include/TriggerFileConfigEntryReader.h"

namespace trigger_tools {

TriggerFileConfigEntryReader::TriggerFileConfigEntryReader(TriggerFileDescriptorCollection& _descriptors)
    : descriptors(&_descriptors) {}

void TriggerFileConfigEntryReader::StartEntry(const std::string& name, const std::string& reference_name)
{
    ConfigEntryReader::StartEntry(name, reference_name);
    current = reference_name.size() ? descriptors->at(reference_name) : TriggerFileDescriptor();
    current.name = name;
}

void TriggerFileConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("channels", 1, Condition::less_equal);
    CheckReadParamCounts("leg", 0, Condition::greater_equal);

    (*descriptors)[current.name] = current;
}

void TriggerFileConfigEntryReader::ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                                            std::istringstream& /*ss*/)
{
    ParseEnumList("channels", current.channels);
    ParseEntry("leg", current.legs);
}

} // namespace trigger_tools
