/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TriggerTools.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "h-tautau/Core/include/TriggerFileDescriptor.h"
#include "h-tautau/Core/include/TriggerFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"

namespace trigger_tools {

    inline constexpr int GetCMSSWVersion()
    {
        int d1 =  *(PROJECT_VERSION + 6) - '0';
        int d2 =  *(PROJECT_VERSION + 7) - '0';
        if(d2 >= 0 && d2 <= 9) return d1 * 10 + d2;
        return d1;
    }

    namespace detail {
        template<typename TriggerObject, int cmssw_version>
        struct UnpackFiltersImpl;

        template<typename TriggerObject>
        struct UnpackFiltersImpl<TriggerObject, 8>{
            static void Unpack(const edm::Event& event, const edm::TriggerResults& triggerResults,
                               TriggerObject& triggerObject) {}
        };

        template<typename TriggerObject>
        struct UnpackFiltersImpl<TriggerObject, 9>{
            static void Unpack(const edm::Event& event, const edm::TriggerResults& triggerResults,
                               TriggerObject& triggerObject)
            {
                triggerObject.unpackFilterLabels(event, triggerResults);
            }
        };

        template<typename TriggerObject>
        struct UnpackFiltersImpl<TriggerObject, 10>{
            static void Unpack(const edm::Event& event, const edm::TriggerResults& triggerResults,
                               TriggerObject& triggerObject)
            {
                triggerObject.unpackFilterLabels(event, triggerResults);
            }
        };

        template<typename TriggerObject>
        struct UnpackFiltersImpl<TriggerObject, 11>{
            static void Unpack(const edm::Event& event, const edm::TriggerResults& triggerResults,
                               TriggerObject& triggerObject)
            {
                triggerObject.unpackFilterLabels(event, triggerResults);
            }
        };
    }

    inline void UnpackFilters(const edm::Event& event, const edm::TriggerResults& triggerResults,
                              pat::TriggerObjectStandAlone& triggerObject)
    {
        detail::UnpackFiltersImpl<pat::TriggerObjectStandAlone, GetCMSSWVersion()>::Unpack(event, triggerResults, triggerObject);
    }
}

namespace analysis {

TriggerTools::TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIMembedding_token,
                           EDGetTokenT<edm::TriggerResults>&& _triggerResultsMERGE_token,
                           EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                           EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
                           EDGetTokenT<BXVector<l1t::Tau>>&& _l1Tau_token,
                           const std::string& triggerCfg, Channel _channel, bool _isEmbedded) :
    triggerPrescales_token(_triggerPrescales_token), triggerObjects_token(_triggerObjects_token),
    l1Tau_token(_l1Tau_token), isEmbedded(_isEmbedded), channel(_channel)
{
    triggerResults_tokens[CMSSW_Process::SIM] = _triggerResultsSIM_token;
    triggerResults_tokens[CMSSW_Process::HLT] = _triggerResultsHLT_token;
    triggerResults_tokens[CMSSW_Process::RECO] = _triggerResultsRECO_token;
    triggerResults_tokens[CMSSW_Process::PAT] = _triggerResultsPAT_token;
    triggerResults_tokens[CMSSW_Process::SIMembedding] = _triggerResultsSIMembedding_token;
    triggerResults_tokens[CMSSW_Process::MERGE] = _triggerResultsMERGE_token;

    triggerDescriptors = TriggerDescriptorCollection::Load(triggerCfg,channel);

    pathTriggerObjects.resize(triggerDescriptors->size());
    for (size_t n = 0; n < pathTriggerObjects.size(); ++n) {
        auto& legId_triggerObjPtr_vector = pathTriggerObjects.at(n);
        legId_triggerObjPtr_vector.resize(triggerDescriptors->at(n).lepton_legs.size(), {});
    }
}

void TriggerTools::Initialize(const edm::Event &_iEvent, bool isData)
{
    iEvent = &_iEvent;
    for(const auto& entry : triggerResults_tokens)
        iEvent->getByToken(entry.second, triggerResultsMap[entry.first]);
    iEvent->getByToken(triggerPrescales_token, triggerPrescales);
    iEvent->getByToken(triggerObjects_token, triggerObjects);
    iEvent->getByToken(l1Tau_token, l1Taus);

    for (auto& desc : pathTriggerObjects){
        for(auto& leg : desc){
            leg.clear();
        }
    }
    jetTriggerObjects = TriggerDescriptorCollection::JetTriggerObjectCollection();

    static const std::map<LegType, std::set<trigger::TriggerObjectType>> map_legType_triggerObjType = {
        { LegType::e, { trigger::TriggerElectron, trigger::TriggerCluster } },
        { LegType::mu, { trigger::TriggerMuon } },
        { LegType::tau, { trigger::TriggerTau } },
        { LegType::jet, { trigger::TriggerJet } }
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

    Handle<edm::TriggerResults>& triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::HLT);
    if(isEmbedded) triggerResultsHLT = triggerResultsMap.at(CMSSW_Process::SIMembedding);

    const edm::TriggerNames& triggerNames = iEvent->triggerNames(*triggerResultsHLT);
    for (const pat::TriggerObjectStandAlone& triggerObject : *triggerObjects) {
        pat::TriggerObjectStandAlone unpackedTriggerObject(triggerObject);
        unpackedTriggerObject.unpackPathNames(triggerNames);
        trigger_tools::UnpackFilters(*iEvent,*triggerResultsHLT,unpackedTriggerObject);
        const auto& paths = unpackedTriggerObject.pathNames(false, true); //isLastFilter=false
        if(hasExpectedType(LegType::jet, unpackedTriggerObject)) {
            const auto& jet_filters = triggerDescriptors->GetJetFilters();
            bool jet_added = false;
            size_t jet_index;
            for(size_t filter_index = 0; filter_index < jet_filters.size(); ++filter_index) {
                if(!unpackedTriggerObject.hasFilterLabel(jet_filters.at(filter_index))) continue;
                if(!jet_added) {
                    jet_index = jetTriggerObjects.Add(unpackedTriggerObject.p4());
                    jet_added = true;
                }
                jetTriggerObjects.SetJetFilterMatchBit(filter_index, jet_index, true);
            }
        } else {
            if(debug)
                std::cout << "Trigger object " << LorentzVectorToString(triggerObject.polarP4()) << '\n';
            for(const auto& path : paths) {
                size_t index;
                const bool pattern_match_found = triggerDescriptors->FindPatternMatch(path,index);
                if(debug)
                    std::cout << '\t' << path << " found: " << pattern_match_found << '\n';
                if(!pattern_match_found) continue;
                const auto& descriptor = triggerDescriptors->at(index);
                for(unsigned n = 0; n < descriptor.lepton_legs.size(); ++n){
                    const TriggerDescriptorCollection::Leg& leg = descriptor.lepton_legs.at(n);
                    TriggerDescriptorCollection::FilterVector filter_toPass;
                    if(isData && leg.run_switch.is_initialized() && _iEvent.id().run() < *leg.run_switch)
                        filter_toPass = *leg.legacy_filters;
                    else
                        filter_toPass = leg.filters;
                    const bool has_expected_type = hasExpectedType(leg.type, unpackedTriggerObject);
                    const bool pass_filters = passFilters(unpackedTriggerObject,filter_toPass);
                    if(debug) {
                        std::cout << "\t\tleg=" << n << ", has_expected_type=" << has_expected_type
                                  << ", pass_filters=" << pass_filters << '\n';
                        for(const auto& filter : filter_toPass) {
                            std::cout << "\t\t\t" << filter << ": "
                                      << unpackedTriggerObject.hasFilterLabel(filter) << '\n';
                        }
                    }
                    if(!has_expected_type || !pass_filters) continue;
                    pathTriggerObjects.at(index).at(n).insert(&triggerObject);
                }
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
        if(triggerDescriptors->FindPatternMatch(triggerNames.triggerName(i), index))
            results.SetAccept(index, triggerResultsHLT->accept(i));
    }
}

TriggerTools::VectorTriggerObjectSet TriggerTools::FindMatchingTriggerObjects(
        size_t index, const LorentzVector& candidateMomentum, LegType candidate_type, double deltaR_Limit) const
{
    const auto& legId_triggerObjPtr_vector = pathTriggerObjects.at(index);
    TriggerTools::VectorTriggerObjectSet matched_legId_triggerObjectSet_vector(legId_triggerObjPtr_vector.size());
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    const auto& descriptor = triggerDescriptors->at(index);

    for(size_t n = 0; n < legId_triggerObjPtr_vector.size(); ++n) {
        const auto& triggerObjectSet = legId_triggerObjPtr_vector.at(n);
        const TriggerDescriptorCollection::Leg& leg = descriptor.lepton_legs.at(n);
        if(candidate_type != leg.type) continue;
        for(const auto& triggerObject : triggerObjectSet) {
            if(ROOT::Math::VectorUtil::DeltaR2(triggerObject->polarP4(), candidateMomentum) >= deltaR2) continue;
            if(leg.eta.is_initialized() && std::abs(candidateMomentum.Eta()) >= leg.eta ) continue;
            //if(candidateMomentum.Pt() <= leg.pt + deltaPt_map.at(leg.type)) continue;
            static const double deltaR2_l1 = std::pow(0.5, 2);
            bool found_l1_match = false;
            if(leg.applyL1match){
                const BXVector<l1t::Tau>& l1taus_elements = *l1Taus.product();
                for (unsigned n = 0; n < l1taus_elements.size(0) && !found_l1_match; ++n){
                    const l1t::Tau& l1tau = l1taus_elements.at(0,n);
                    found_l1_match = l1tau.hwIso() > 0.5 && l1tau.et() > 32 && leg.type == analysis::LegType::tau &&
                                     ROOT::Math::VectorUtil::DeltaR2(l1tau.p4(), candidateMomentum) < deltaR2_l1;
                }
            }
            if (!leg.applyL1match || found_l1_match) matched_legId_triggerObjectSet_vector.at(n).insert(triggerObject);
        }
    }

    return matched_legId_triggerObjectSet_vector;
}

bool TriggerTools::TriggerMatchFound(const std::array<TriggerTools::VectorTriggerObjectSet, 2>& matched_legIds,
                       const size_t n_legs_total)
{
    if(n_legs_total == 0) return true;

    if(n_legs_total == 1)
        return matched_legIds.at(0).at(0).size() >= n_legs_total || matched_legIds.at(1).at(0).size() >= n_legs_total;

    bool match_found = false;
    if(n_legs_total == 2) {
        for(size_t flip = 0; !match_found && flip < matched_legIds.size(); ++flip) {
            const size_t first = flip, second = ((flip + 1) % 2);
            std::vector<const pat::TriggerObjectStandAlone*> comb_match;
            std::set_union(matched_legIds.at(0).at(first).begin(),
                           matched_legIds.at(0).at(first).end(),
                           matched_legIds.at(1).at(second).begin(),
                           matched_legIds.at(1).at(second).end(),
                           std::back_inserter(comb_match));

            match_found = matched_legIds.at(0).at(first).size() >= 1 &&
                    matched_legIds.at(1).at(second).size() >= n_legs_total - 1 &&
                    comb_match.size() >= n_legs_total;
        }
    }
    return match_found;
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
    std::map<CMSSW_Process, bool> results;
    for(const auto& process : all_processes) {
        if(TryGetTriggerResult(process, name, result)) {
            if(!debug)
                return true;
            results[process] = result;
        }
    }
    if(!debug || results.empty())
        return false;
    if(results.size() > 1) {
        for(const auto& entry : results)
            std::cout << entry.first << ": " << entry.second << std::endl;
        throw analysis::exception("Multiple trigger results for '%1%'.") % name;
    }
    return true;
}

bool TriggerTools::GetAnyTriggerResult(const std::string& name) const
{
    bool result = false;
    if(TryGetAnyTriggerResult(name, result))
        return result;
    throw analysis::exception("Unable to find trigger '%1%'") % name;
}

TriggerTools::BitsContainer TriggerTools::GetJetMatchBitsImpl(const LorentzVector& reco_jet_p4,
                                                              double deltaR_Limit) const
{
    auto trig_objs = jetTriggerObjects;
    const auto& jet_filters = triggerDescriptors->GetJetFilters();
    const double deltaR2 = std::pow(deltaR_Limit, 2);
    for(size_t filter_index = 0; filter_index < jet_filters.size(); ++filter_index) {
        for(size_t jet_index = 0; jet_index < trig_objs.momentums.size(); ++jet_index) {
            const bool match = trig_objs.GetJetFilterMatchBit(filter_index, jet_index)
                && ROOT::Math::VectorUtil::DeltaR2(trig_objs.momentums.at(jet_index), reco_jet_p4) < deltaR2;
            trig_objs.SetJetFilterMatchBit(filter_index, jet_index, match);
        }
    }
    return trig_objs.match_bits;
}


} // namespace analysis
