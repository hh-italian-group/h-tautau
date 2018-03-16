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
                           std::string triggerCfg, Channel channel) :
    triggerPrescales_token(_triggerPrescales_token), triggerObjects_token(_triggerObjects_token),
    l1JetParticles_token(_l1JetParticles_token)
{
    triggerResults_tokens[CMSSW_Process::SIM] = _triggerResultsSIM_token;
    triggerResults_tokens[CMSSW_Process::HLT] = _triggerResultsHLT_token;
    triggerResults_tokens[CMSSW_Process::RECO] = _triggerResultsRECO_token;
    triggerResults_tokens[CMSSW_Process::PAT] = _triggerResultsPAT_token;

    trigger_tools::SetupDescriptor setup;
    trigger_tools::TriggerFileDescriptorCollection trigger_file_descriptors = TriggerTools::ReadConfig(triggerCfg,setup);

    deltaPt_map = setup.deltaPt_map;

    triggerDescriptors = CreateTriggerDescriptors(trigger_file_descriptors,channel);

    path_legId_triggerObjPtr_vector.resize(triggerDescriptors.size());
    for (size_t n = 0; n < path_legId_triggerObjPtr_vector.size(); ++n){
        const auto& legId_triggerObjPtr_vector = path_legId_triggerObjPtr_vector.at(n);
        legId_triggerObjPtr_vector.resize(triggerDescriptors.at(n).legs_info.size(), {});
    }

}

trigger_tools::TriggerFileDescriptorCollection ReadConfig(const std::string& cfg_path,
                                                          trigger_tools::SetupDescriptor& setup)
{
    trigger_tools::TriggerFileDescriptorCollection trigger_file_descriptors;
    analysis::ConfigReader config_reader;
    trigger_tools::TriggerFileConfigEntryReader trigger_entry_reader(trigger_file_descriptors);
    config_reader.AddEntryReader("PATTERN", trigger_entry_reader, true);

    trigger_tools::SetupDescriptorCollection setup_file_descriptors;
    trigger_tools::SetupConfigEntryReader setup_entry_reader(setup_file_descriptors);
    config_reader.AddEntryReader("SETUP", setup_entry_reader, false);

    const std::string triggerCfg_full = edm::FileInPath(cfg_path).fullPath();
    config_reader.ReadConfig(triggerCfg_full);

    if(setup_file_descriptors.size() != 1)
        throw exception("More than 1 setup in Reading Trigger Tools cfg");
    setup = *setup_file_descriptors.begin();

    return trigger_file_descriptors;
}

TriggerDescriptorCollection TriggerTools::CreateTriggerDescriptors(const trigger_tools::TriggerFileDescriptorCollection& trigger_file_descriptors,
                                                                         Channel channel)
{
    TriggerDescriptorCollection triggerDescriptors;
    for(const auto& entry : trigger_file_descriptors) {
        trigger_tools::TriggerFileDescriptor trigger_file_descriptor = entry.second;
        if(!trigger_file_descriptor.channels.count(channel)) continue;
        const auto& legs = trigger_file_descriptor.legs;
        std::vector<TriggerDescriptorCollection::Leg> legs_vector;
        for (size_t n = 0; n < legs.size(); ++n){
            const analysis::PropertyList leg_list = analysis::Parse<analysis::PropertyList>(legs.at(n));
            const analysis::LegType type = leg_list.Get<analysis::LegType>("type");
            const double pt = leg_list.Get<double>("pt");
            const TriggerDescriptorCollection::FilterVector filters = leg_list.GetList<std::string>("filters", false);
            legs_vector.emplace_back(type,pt,filters);
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

    for (auto& desc : path_legId_triggerObjPtr_vector){
        for(auto& leg : desc){
            leg.clear();
        }
    }

    static const std::map<LegType,std::set<trigger::TriggerObjectType>> map_legType_triggerObjType
            = {
               { LegType::e, { trigger::TriggerElectron, trigger::TriggerCluster } },
               { LegType::mu, { trigger::TriggerMuon } },
               { LegType::tau, { trigger::TriggerTau } }
              };

    const auto hasExpectedType = [](LegType legType, const pat::TriggerObjectStandAlone& triggerObject) {
            for(auto type : map_legType_triggerObjType.at(legType))
                if(triggerObject.type(type)) return true;
            return false;
    };


    const auto passFilters = [](const pat::TriggerObjectStandAlone& triggerObject,
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
            size_t index;
            if(!triggerDescriptors.FindPatternMatch(path,index)) continue;
            const auto& descriptor = triggerDescriptors.at(index);
            for(unsigned n = 0; n < descriptor.legs_info.size(); ++n){
                const TriggerDescriptorCollection::Leg& leg = descriptor.legs_info.at(n);
                if(!hasExpectedType(leg.type,unpackedTriggerObject) || !passFilters(unpackedTriggerObject,leg.filters)) continue;
                path_legId_triggerObjPtr_vector.at(index).at(n).insert(&triggerObject);
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

std::vector<TriggerTools::TriggerObjectSet> FindMatchingTriggerObjects(
        const size_t index, const LorentzVector& candidateMomentum, const LegType& candidate_type, double deltaR_Limit)
{
    const auto& legId_triggerObjPtr_vector = path_legId_triggerObjPtr_vector.at(index);
    std::vector<TriggerTools::TriggerObjectSet> matched_legId_triggerObjectSet_vector(legId_triggerObjPtr_vector.size());
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    
    for(size_t n= 0; n < legId_triggerObjPtr_vector.size(); ++n){
        const auto& triggerObjectSet = legId_triggerObjPtr_vector.at(n);
        const auto& descriptor = triggerDescriptors.at(index);
        const TriggerDescriptorCollection::Leg& leg = descriptor.legs_info.at(n);
        if(candidate_type != leg.type) continue;
        for(const auto& triggerObject : triggerObjectSet){
            if(ROOT::Math::VectorUtil::DeltaR2(triggerObject->polarP4(), candidateMomentum) >= deltaR2) continue;
            if(candidateMomentum.Pt() <= leg.pt + deltaPt_map.at(leg.type)) continue;
            matched_legId_triggerObjectSet_vector.at(n).insert(triggerObject);
        }
    }

    return matched_legId_triggerObjectSet_vector;
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
