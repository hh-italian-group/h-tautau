/*! Definition of the file configuration entry reader for DY sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TriggerFileDescriptor.h"

namespace analysis {

namespace trigger{

class TriggerFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TriggerFileConfigEntryReader(TriggerDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : NJets_HT_BinFileDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("n_jet", 1, Condition::equal_to);
        CheckReadParamCounts("n_bjet", 1, Condition::equal_to);
        CheckReadParamCounts("n_ht", 1, Condition::equal_to);
        CheckReadParamCounts("fileType", 1, Condition::equal_to);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("n_jet", current.n_jet);
        ParseEntry("n_bjet", current.n_bjet);
        ParseEntry("n_ht", current.n_ht);
        ParseEntry("fileType", current.fileType);
    }

private:
    TriggerFileDescriptor current;
    TriggerDescriptorCollection* descriptors;
};

} //namespace trigger

} // namespace analysis
