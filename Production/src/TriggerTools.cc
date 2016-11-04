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
}

bool TriggerTools::HaveTriggerFired(const std::vector<std::string>& hltPaths)
{
    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);

    for (unsigned i = 0; i < triggerResultsHLT->size(); ++i) {
        if(triggerPrescales->getPrescaleForIndex(i) != 1 && triggerResultsHLT->accept(i)) continue;
        const std::string& objectMatchedPath = triggerNames.triggerName(i);
        for (const std::string& triggerPath : hltPaths ){
            const size_t found = objectMatchedPath.find(triggerPath);
            if(found != std::string::npos)
                return true;
        }
    }
    return false;
}

std::set<const pat::TriggerObjectStandAlone*> TriggerTools::FindMatchingTriggerObjects(
        const std::string& pathOfInterest, trigger::TriggerObjectType objectType,
        const LorentzVector& candidateMomentum, double deltaR_Limit)
{
    std::set<const pat::TriggerObjectStandAlone*> matches;
    const auto& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        if(!triggerObject.type(objectType)) continue;
        if(ROOT::Math::VectorUtil::DeltaR2(triggerObject.polarP4(), candidateMomentum) >= deltaR2) continue;
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        const auto& matchedPaths = unpackedTriggerObject.pathNames(true, true);
        for(const auto& matchedPath : matchedPaths) {
            if( matchedPath.find(pathOfInterest) == std::string::npos ) continue;
//            const bool LF = unpackedTriggerObject.hasPathName(matchedPath, true, false);
//            const bool L3 = unpackedTriggerObject.hasPathName(matchedPath, false, true);
//            const bool all = unpackedTriggerObject.hasPathName(matchedPath, true, true);
//            std::cout << "RECO p4 = " << ConvertVector(candidateMomentum)
//                      << ", trig obj p4 = " << ConvertVector(triggerObject.polarP4())
//                      << ", deltaR = " << ROOT::Math::VectorUtil::DeltaR(triggerObject.polarP4(), candidateMomentum)
//                      << ", matchedPath = " << matchedPath
//                      << ", LF = " << LF
//                      << ", L3 = " << L3
//                      << ", all = " << all
//                      << ", obj types = ( ";
//            for(int type : triggerObject.triggerObjectTypes())
//                std::cout << type << " ";
//            std::cout << ")" << std::endl;

            matches.insert(&triggerObject);
            break;
        }
    }
    return matches;
}

TriggerTools::L1ParticlePtrSet TriggerTools::L1TauMatch(const LorentzVector& tauMomentum)
{
    static constexpr double MinPt = 28;
    static constexpr double DeltaR = 0.5;
    static constexpr double DeltaR2 = std::pow(DeltaR, 2);

    L1ParticlePtrSet result;
    for(const auto& l1Tau : *l1JetParticles){
        if(l1Tau.pt() > MinPt && ROOT::Math::VectorUtil::DeltaR2(l1Tau.polarP4(), tauMomentum) < DeltaR2)
            result.insert(&l1Tau);
    }
    return result;
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

bool TriggerTools::GetAnyTriggerResult(const std::string& name) const
{
    static const auto all_processes = __CMSSW_Process_names<>::names.GetEnumEntries();
    for(const auto& process : all_processes) {
        bool result;
        if(TryGetTriggerResult(process, name, result))
            return result;
    }
    throw analysis::exception("Unable to find trigger '%1%'") % name;
}

} // namespace analysis
