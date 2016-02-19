/*! Definiton of tools to work with embedded trigger information.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/PatCandidates/interface/PATObject.h"

template<typename PatObject>
std::vector<std::string> CollectMatchedTriggerPaths(const PatObject& patObject)
{
    const pat::TriggerObjectStandAloneCollection& matchedTriggers = patObject.triggerObjectMatches();
    std::set<std::string> pathNames;
    for(const pat::TriggerObjectStandAlone& triggerObject : matchedTriggers) {
        const auto objectPathNames = triggerObject.pathNames();
        for(const std::string& pathName : objectPathNames)
            pathNames.insert(pathName);
    }
    return std::vector<std::string>(pathNames.begin(), pathNames.end());
}
