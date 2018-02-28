/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TriggerTools.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {

TriggerTools::TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                           EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                           EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                           EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token) :
    triggerPrescales_token(_triggerPrescales_token), triggerObjects_token(_triggerObjects_token),
    l1JetParticles_token(_l1JetParticles_token)
{
    triggerResults_tokens[CMSSW_Process::SIM] = _triggerResultsSIM_token;
    triggerResults_tokens[CMSSW_Process::HLT] = _triggerResultsHLT_token;
    triggerResults_tokens[CMSSW_Process::RECO] = _triggerResultsRECO_token;
    triggerResults_tokens[CMSSW_Process::PAT] = _triggerResultsPAT_token;
}

void TriggerTools::Initialize(const edm::Event &_iEvent)
{
    iEvent = &_iEvent;
    for(const auto& entry : triggerResults_tokens)
        iEvent->getByToken(entry.second, triggerResultsMap[entry.first]);
    iEvent->getByToken(triggerPrescales_token, triggerPrescales);
    iEvent->getByToken(triggerObjects_token, triggerObjects);
    iEvent->getByToken(l1JetParticles_token, l1JetParticles);

    const auto hasExpectedType = [&](const pat::TriggerObjectStandAlone& triggerObject) {
        for(auto type : objectTypes)
            if(triggerObject.type(type)) return true;
        return false;
    };

    const auto& filters = descriptors.GetFilters(path_index, leg_id);
    const auto passFilters = [&](const pat::TriggerObjectStandAlone& triggerObject) {
        for(const auto& filter : filters)
            if(!triggerObject.hasFilterLabel(filter)) return false;
        return true;
    };

    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);

    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        if(!hasExpectedType(triggerObject)) continue;
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        unpackedTriggerObject.unpackFilterLabels(*iEvent,*triggerResultsHLT); //new
        if(!passFilters(unpackedTriggerObject)) continue;
        const auto& paths = unpackedTriggerObject.pathNames(true, true);
        for(const auto& path : paths) {
            path_triggerObjPtr_map[path].push_back(*triggerObject);
        }
    }

}

void TriggerTools::SetTriggerAcceptBits(const analysis::TriggerDescriptors& descriptors,
                                        analysis::TriggerResults& results)
{
    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);

    for (size_t i = 0; i < triggerResultsHLT->size(); ++i) {
        if(triggerPrescales->getPrescaleForIndex(i) != 1) continue;
        size_t index;
        if(descriptors.FindPatternMatch(triggerNames.triggerName(i), index))
            results.SetAccept(index, triggerResultsHLT->accept(i));
    }
}

TriggerTools::TriggerObjectSet TriggerTools::FindMatchingTriggerObjects(
        const TriggerDescriptors& descriptors, size_t path_index,
        const std::set<trigger::TriggerObjectType>& objectTypes, const LorentzVector& candidateMomentum,
        size_t leg_id, double deltaR_Limit)
{
    const auto hasExpectedType = [&](const pat::TriggerObjectStandAlone& triggerObject) {
        for(auto type : objectTypes)
            if(triggerObject.type(type)) return true;
        return false;
    };

    const auto& filters = descriptors.GetFilters(path_index, leg_id);
    const auto passFilters = [&](const pat::TriggerObjectStandAlone& triggerObject) {
        for(const auto& filter : filters)
            if(!triggerObject.hasFilterLabel(filter)) return false;
        return true;
    };
    
    TriggerObjectSet matches;
    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        if(!hasExpectedType(triggerObject)) continue;
        if(ROOT::Math::VectorUtil::DeltaR2(triggerObject.polarP4(), candidateMomentum) >= deltaR2) continue;
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        unpackedTriggerObject.unpackFilterLabels(*iEvent,*triggerResultsHLT); //new
        if(!passFilters(unpackedTriggerObject)) continue;
        const auto& paths = unpackedTriggerObject.pathNames(true, true);
        for(const auto& path : paths) {
            if(descriptors.PatternMatch(path, path_index)) {
                matches.insert(&triggerObject);
                break;
            }
        }
    }
    return matches;
}

bool TriggerTools::TryGetTriggerResult(CMSSW_Process process, const std::string& name, bool& result) const
{
    const auto& triggerResults = triggerResultsMap.at(process);
    if(!triggerResults.isValid()) return false;
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResults);
    size_t index = triggerNames.triggerIndex(name);
    if(index == triggerNames.size()) return false;
    result = triggerResults->accept(index);
    return true;
}

bool TriggerTools::GetTriggerResult(CMSSW_Process process, const std::string& name) const
{
    bool result;
    if(!TryGetTriggerResult(process, name, result)) {
        std::ostringstream ss;
        ss << "Unable to find trigger '" << name << "'. ";
        const auto& triggerResults = triggerResultsMap.at(process);
        if(!triggerResults.isValid())
            ss << "Input collection for " << process << " process not found.";
        else {
            ss << "Available triggers:\n";
            const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResults);
            for(const auto& name : triggerNames.triggerNames())
                ss << "\t" << name << "\n";
        }
        throw analysis::exception(ss.str());
    }
    return result;
}

bool TriggerTools::TryGetAnyTriggerResult(const std::string& name, bool& result) const
{
    static const auto& all_processes = __CMSSW_Process_names<>::names.GetEnumEntries();
    for(auto process : all_processes) {
        if(TryGetTriggerResult(process, name, result))
            return true;
    }
    return false;
}

bool TriggerTools::GetAnyTriggerResult(const std::string& name) const
{
    bool result;
    if(TryGetAnyTriggerResult(name, result))
        return result;
    throw analysis::exception("Unable to find trigger '%1%'") % name;
}

} // namespace analysis
