/*! Definition of the file configuration entry reader for Trigger configuration file.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TriggerFileDescriptor.h"

namespace trigger_tools {

class TriggerFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TriggerFileConfigEntryReader(TriggerFileDescriptorCollection& _descriptors);

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override;
    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;

private:
    TriggerFileDescriptor current;
    TriggerFileDescriptorCollection* descriptors;
};

class SetupConfigEntryReader : public analysis::ConfigEntryReaderT<SetupDescriptor> {
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
