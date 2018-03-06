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
                           EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token,
                           std::vector<edm::ParameterSet> _hltPaths) :
    triggerPrescales_token(_triggerPrescales_token), triggerObjects_token(_triggerObjects_token),
    l1JetParticles_token(_l1JetParticles_token), hltPaths(_hltPaths)
{
    triggerResults_tokens[CMSSW_Process::SIM] = _triggerResultsSIM_token;
    triggerResults_tokens[CMSSW_Process::HLT] = _triggerResultsHLT_token;
    triggerResults_tokens[CMSSW_Process::RECO] = _triggerResultsRECO_token;
    triggerResults_tokens[CMSSW_Process::PAT] = _triggerResultsPAT_token;

    for(const auto& hltPath : hltPaths) {
        pattern = hltPath.getParameter<std::string>("pattern");
        legs = hltPath.getUntrackedParameter<std::vector<std::string>>("legs", {});
        TriggerTools::CreateTriggerDescriptors();
    }

}

void TriggerTools::CreateTriggerDescriptors()
{
    std::vector<analysis::TriggerDescriptors::Leg> legs_vector;
    for (unsigned n = 0; n < legs.size(); ++n){
        const analysis::PropertyList leg_list = analysis::Parse<analysis::PropertyList>(legs.at(n));
        const analysis::LegType type = leg_list.Get<analysis::LegType>("type");
        const double pt = leg_list.Get<double>("pt");
        const analysis::TriggerDescriptors::FilterVector filters = analysis::PropertyList::GetList(leg_list.Get("filters"),false);
        const analysis::TriggerDescriptors::Leg legs_struct(type,pt,filters);
        legs_vector.push_back(legs_struct);
    }
    triggerDescriptors.Add(pattern, legs_vector);
}

void TriggerTools::Initialize(const edm::Event &_iEvent)
{
    iEvent = &_iEvent;
    for(const auto& entry : triggerResults_tokens)
        iEvent->getByToken(entry.second, triggerResultsMap[entry.first]);
    iEvent->getByToken(triggerPrescales_token, triggerPrescales);
    iEvent->getByToken(triggerObjects_token, triggerObjects);
    iEvent->getByToken(l1JetParticles_token, l1JetParticles);

    const auto hasExpectedType = [&](const LegType& legType, const pat::TriggerObjectStandAlone& triggerObject) {
            for(auto type : detail::GetTriggerObjectTypes().at(legType))
                if(triggerObject.type(type)) return true;
            return false;
    };


    const auto passFilters = [&](const pat::TriggerObjectStandAlone& triggerObject,
            const TriggerDescriptors::FilterVector& filters) {
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
            const TriggerDescriptors::PatternStruct pattern_struct = triggerDescriptors.GetPatternStruct(path);
            for(unsigned n = 0; n < pattern_struct.legs_info.size(); ++n){
                const TriggerDescriptors::Leg leg = pattern_struct.legs_info.at(n);
                if(!passFilters(unpackedTriggerObject,leg.filters)) continue;
                if(!hasExpectedType(leg.type,unpackedTriggerObject)) continue;
                path_legId_triggerObjPtr_map[pattern_struct.pattern][n].insert(*unpackedTriggerObject);
            }
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
        const TriggerDescriptors& descriptors, const LorentzVector& candidateMomentum,
        double deltaR_Limit)
{



    
    TriggerObjectSet matches;

    const double deltaR2 = std::pow(deltaR_Limit, 2);

    //loop on map
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {

        if(ROOT::Math::VectorUtil::DeltaR2(triggerObject.polarP4(), candidateMomentum) >= deltaR2) continue;


            if(descriptors.PatternMatch(path, path_index)) {
                matches.insert(&triggerObject);
                break;
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
