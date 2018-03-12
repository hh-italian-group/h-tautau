/*! Definition of the file configuration entry reader for DY sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TriggerFileDescriptor.h"

namespace analysis {

namespace trigger{

class TriggerFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TriggerFileConfigEntryReader(TriggerFileDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : TriggerFileDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("leg", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEnumList("channels", current.channels);
        ParseEntry("leg", current.legs);

    }

private:
    analysis::trigger::TriggerFileDescriptor current;
    analysis::trigger::TriggerFileDescriptorCollection* descriptors;
};

class SetupConfigEntryReader : public analysis::ConfigEntryReader {
public:
    SetupConfigEntryReader(SetupDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : SetupDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("deltaPt", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("deltaPt", current.deltaPt_map);
    }

private:
    analysis::trigger::SetupDescriptor current;
    analysis::trigger::SetupDescriptorCollection* descriptors;
};

} //namespace trigger

} // namespace analysis
