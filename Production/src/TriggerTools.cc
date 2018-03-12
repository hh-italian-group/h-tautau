/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TriggerTools.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "h-tautau/Production/interface/TriggerFileDescriptor.h"
#include "h-tautau/Production/interface/TriggerFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"

namespace analysis {

TriggerTools::TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                           EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                           EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                           EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token,
                           std::string _triggerCfg,
                           analysis::Channel _channel) :
    triggerPrescales_token(_triggerPrescales_token), triggerObjects_token(_triggerObjects_token),
    l1JetParticles_token(_l1JetParticles_token), triggerCfg(_triggerCfg), channel(_channel)
{
    triggerResults_tokens[CMSSW_Process::SIM] = _triggerResultsSIM_token;
    triggerResults_tokens[CMSSW_Process::HLT] = _triggerResultsHLT_token;
    triggerResults_tokens[CMSSW_Process::RECO] = _triggerResultsRECO_token;
    triggerResults_tokens[CMSSW_Process::PAT] = _triggerResultsPAT_token;

    ConfigReader config_reader;

    TriggerFileDescriptorCollection trigger_file_descriptors;
    TriggerFileConfigEntryReader trigger_entry_reader(trigger_file_descriptors);
    config_reader.AddEntryReader("PATTERN", trigger_entry_reader, true);

    SetupDescriptorCollection setup_descriptors;
    SetupConfigEntryReader setup_entry_reader(setup_descriptors);
    config_reader.AddEntryReader("SETUP", setup_entry_reader, false);

    config_reader.ReadConfig(triggerCfg);

    for (const auto& setup : setup_descriptors){
        const SetupDescriptor setup_descriptor = setup.second;
        for(const auto iter : setup_descriptor.deltaPt_map){
            deltaPt_map[iter.first] = iter.second;
        }
    }

    triggerDescriptors = TriggerTools::CreateTriggerDescriptors(trigger_file_descriptors,channel);

}

TriggerDescriptorCollection TriggerTools::CreateTriggerDescriptors(const analysis::TriggerFileDescriptorCollection trigger_file_descriptors,
                                                                          const Channel channel)
{
    TriggerDescriptorCollection triggerDescriptors;
    for(const auto& entry : trigger_file_descriptors) {
        TriggerFileDescriptor trigger_file_descriptor = entry.second;
        if(!trigger_file_descriptor.channels.count(channel)) continue;
        const std::vector<std::string> legs = trigger_file_descriptor.legs;
        std::vector<analysis::TriggerDescriptorCollection::Leg> legs_vector;
        for (unsigned n = 0; n < legs.size(); ++n){
            const analysis::PropertyList leg_list = analysis::Parse<analysis::PropertyList>(legs.at(n));
            const analysis::LegType type = leg_list.Get<analysis::LegType>("type");
            const double pt = leg_list.Get<double>("pt");
            const analysis::TriggerDescriptorCollection::FilterVector filters = leg_list.GetList(leg_list.Get("filters"),false);
            const analysis::TriggerDescriptorCollection::Leg legs_struct(type,pt,filters);
            legs_vector.push_back(legs_struct);
        }
        triggerDescriptors.Add(entry.first, legs_vector);
    }
    return triggerDescriptors;
}

void TriggerTools::Initialize(const edm::Event &_iEvent)
{
    iEvent = &_iEvent;
    for(const auto& entry : triggerResults_tokens)
        iEvent->getByToken(entry.second, triggerResultsMap[entry.first]);
    iEvent->getByToken(triggerPrescales_token, triggerPrescales);
    iEvent->getByToken(triggerObjects_token, triggerObjects);
    iEvent->getByToken(l1JetParticles_token, l1JetParticles);

    static const std::map<analysis::LegType,std::set<trigger::TriggerObjectType>> map_legType_triggerObjType
            = {
               { analysis::LegType::e, { trigger::TriggerElectron, trigger::TriggerCluster } },
               { analysis::LegType::mu, { trigger::TriggerMuon } },
               { analysis::LegType::tau, { trigger::TriggerTau } }
              };

    const auto hasExpectedType = [&](const LegType& legType, const pat::TriggerObjectStandAlone& triggerObject) {
            for(auto type : map_legType_triggerObjType.at(legType))
                if(triggerObject.type(type)) return true;
            return false;
    };


    const auto passFilters = [&](const pat::TriggerObjectStandAlone& triggerObject,
            const TriggerDescriptorCollection::FilterVector& filters) {
        for(const auto& filter : filters)
            if(!triggerObject.hasFilterLabel(filter)) return false;
        return true;
    };


    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        unpackedTriggerObject.unpackFilterLabels(*iEvent,*triggerResultsHLT); //new
        const auto& paths = unpackedTriggerObject.pathNames(true, true);
        for(const auto& path : paths) {
            const TriggerDescriptorCollection::TriggerDescriptor pattern_struct = triggerDescriptors.GetTriggerDescriptor(path);
            for(unsigned n = 0; n < pattern_struct.legs_info.size(); ++n){
                const TriggerDescriptorCollection::Leg leg = pattern_struct.legs_info.at(n);
                if(!passFilters(unpackedTriggerObject,leg.filters)) continue;
                if(!hasExpectedType(leg.type,unpackedTriggerObject)) continue;
                path_legId_triggerObjPtr_map[pattern_struct.pattern][n].insert(&triggerObject);
            }
        }
    }

}



void TriggerTools::SetTriggerAcceptBits(analysis::TriggerResults& results)
{
    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);

    for (size_t i = 0; i < triggerResultsHLT->size(); ++i) {
        if(triggerPrescales->getPrescaleForIndex(i) != 1) continue;
        size_t index;
        if(triggerDescriptors.FindPatternMatch(triggerNames.triggerName(i), index))
            results.SetAccept(index, triggerResultsHLT->accept(i));
    }
}

std::map<size_t,TriggerTools::TriggerObjectSet> TriggerTools::FindMatchingTriggerObjects(
        const analysis::TriggerDescriptorCollection::TriggerDescriptor& pattern_struct,
        const LorentzVector& candidateMomentum, const LegType& candidate_type, double deltaR_Limit)
{
    std::map<size_t,TriggerTools::TriggerObjectSet> matched_legId_triggerObjectSet_map;
    const analysis::TriggerDescriptorCollection::Pattern pattern = pattern_struct.pattern;
    std::map<size_t,TriggerTools::TriggerObjectSet> legId_triggerObjPtr_map = path_legId_triggerObjPtr_map.at(pattern);
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    
    for(const auto& iter : legId_triggerObjPtr_map){
        TriggerTools::TriggerObjectSet triggerObjectSet = iter.second;
        for(const auto& triggerObject : triggerObjectSet){
            if(ROOT::Math::VectorUtil::DeltaR2(triggerObject->polarP4(), candidateMomentum) >= deltaR2) continue;
            for(size_t n = 0; n < pattern_struct.legs_info.size(); ++n){
                const TriggerDescriptorCollection::Leg leg = pattern_struct.legs_info.at(n);
                if(candidate_type != leg.type) continue;
                if(candidateMomentum.Pt() <= leg.pt + deltaPt_map.at(leg.type)) continue;
                matched_legId_triggerObjectSet_map[iter.first].insert(*triggerObject);
            }
        }
    }  
    return matched_legId_triggerObjectSet_map;
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
