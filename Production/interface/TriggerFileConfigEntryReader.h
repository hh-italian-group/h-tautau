/*! Definition of the file configuration entry reader for Trigger configuration file.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TriggerFileDescriptor.h"

namespace trigger_tools {

class TriggerFileConfigEntryReader : public ConfigEntryReaderT<TriggerFileDescriptor> {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<TriggerFileDescriptor>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("leg", 0, Condition::greater_equal);

        ConfigEntryReaderT<TriggerFileDescriptor>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEnumList("channels", current.channels);
        ParseEntry("leg", current.legs);

    }
};

class SetupConfigEntryReader : public ConfigEntryReaderT<SetupDescriptor> {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<SetupDescriptor>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("deltaPt", 0, Condition::greater_equal);

        ConfigEntryReaderT<SetupDescriptor>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("deltaPt", current.deltaPt_map);
    }
};

} // namespace analysis
