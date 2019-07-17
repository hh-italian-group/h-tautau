/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/TriggerResults.h"
#include <boost/regex.hpp>
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

TriggerDescriptorCollection::JetTriggerObjectCollection::JetTriggerObjectCollection() : match_bits(0) {}

bool TriggerDescriptorCollection::JetTriggerObjectCollection::GetJetFilterMatchBit(size_t filter_index,
                                                                                   size_t jet_index) const
{
    return ((match_bits >> (filter_index * MaxNumberOfTriggerJets + jet_index)) & BitsContainer(1)) != BitsContainer(0);
}

void TriggerDescriptorCollection::JetTriggerObjectCollection::SetJetFilterMatchBit(size_t filter_index,
                                                                                   size_t jet_index,
                                                                                   bool match_result)
{
    const BitsContainer mask = BitsContainer(1) << (filter_index * MaxNumberOfTriggerJets + jet_index);
    match_bits = match_result ? match_bits | mask : match_bits & ~mask;
}

TriggerDescriptorCollection::Leg::Leg(const LegType _type, double _pt, double _delta_pt, boost::optional<double> _eta, bool _applyL1match,
    const FilterVector& _filters)
    : type(_type), pt(_pt), delta_pt(_delta_pt), eta(_eta), applyL1match(_applyL1match), filters(_filters) { }

TriggerDescriptorCollection::TriggerDescriptor::TriggerDescriptor(const Pattern _pattern,
                                                                  const std::vector<Leg>& legs_info) :
    pattern(_pattern)
{
    static const std::string regex_format = "^%1%[0-9]+$";
    const std::string regex_str = boost::str(boost::format(regex_format) % pattern);
    regex = boost::regex(regex_str);
    for(const auto& leg : legs_info) {
        if(leg.type == LegType::jet)
            jet_legs.emplace_back(leg);
        else
            lepton_legs.emplace_back(leg);
    }
}

bool TriggerDescriptorCollection::TriggerDescriptor::PatternMatch(const std::string& path_name) const
{
    return boost::regex_match(path_name, regex);
}

bool TriggerDescriptorCollection::TriggerDescriptor::RequiresJetMatching() const
{
    return !jet_legs.empty();
}

size_t TriggerDescriptorCollection::size() const { return descriptors.size(); }

const TriggerDescriptorCollection::TriggerDescriptor& TriggerDescriptorCollection::at(size_t index) const
{
    if(index >= descriptors.size())
        throw exception("Index out of trigger descriptors range.");
    return descriptors.at(index);
}

const TriggerDescriptorCollection::TriggerDescriptor& TriggerDescriptorCollection::at(const Pattern& pattern) const
{
    size_t index = GetIndex(pattern);
    return descriptors.at(index);
}

void TriggerDescriptorCollection::Add(const Pattern& pattern, const std::vector<Leg>& legs)
{
    if(desc_indices.count(pattern))
        throw exception("Duplicated trigger pattern '%1%'.") % pattern;
    if(descriptors.size() == TriggerResults::MaxNumberOfTriggers)
        throw exception("Maximal number of triggers is exceeded.");
    desc_indices[pattern] = descriptors.size();
    descriptors.emplace_back(pattern, legs);
    path_index_cache.clear();
    for(auto& leg : descriptors.back().jet_legs) {
        leg.jet_filter_indices.clear();
        for(const auto& filter : leg.filters) {
            auto iter = std::find(jet_filters.begin(), jet_filters.end(), filter);
            if(iter == jet_filters.end()) {
                if(jet_filters.size() == MaxNumberOfJetFilters)
                    throw exception("Maximal number of jet trigger filters is exceeded.");
                jet_filters.push_back(filter);
                iter = std::prev(jet_filters.end());
            }
            leg.jet_filter_indices.push_back(static_cast<unsigned>(std::distance(jet_filters.begin(), iter)));
        }
    }
}

bool TriggerDescriptorCollection::FindPatternMatch(const std::string& path_name, size_t& index)
{
    auto iter = path_index_cache.find(path_name);
    if(iter == path_index_cache.end()) {
        size_t counter = 0;
        for(size_t n = 0; n < descriptors.size(); ++n) {
            const TriggerDescriptor& descriptor = descriptors.at(n);
            if(descriptor.PatternMatch(path_name)) {
                ++counter;
                index = n;
            }
        }
        if (counter > 1)
            throw exception("More than 1 pattern matched.");
        path_index_cache[path_name] = counter == 1 ? index : descriptors.size();
        iter = path_index_cache.find(path_name);
    }
    index = iter->second;
    return index < descriptors.size();
}

size_t TriggerDescriptorCollection::GetIndex(const Pattern& pattern) const
{
    if(!desc_indices.count(pattern))
        throw exception("Trigger '%1%' not found.") % pattern;
    return desc_indices.at(pattern);
}

TriggerDescriptorCollection::FilterBitsContainer TriggerDescriptorCollection::GetJetFilterMatchBits(
    BitsContainer match_bits, unsigned filter_index)
{
    static constexpr FilterBitsContainer mask((FilterBitsContainer(1) << MaxNumberOfTriggerJets) - 1);
    const auto filter_bits = (match_bits >> (filter_index * MaxNumberOfTriggerJets)).convert_to<FilterBitsContainer>();
    return filter_bits & mask;
}

TriggerDescriptorCollection::RootBitsContainer TriggerDescriptorCollection::ConvertToRootRepresentation(
    BitsContainer match_bits)
{
    TriggerDescriptorCollection::RootBitsContainer result;
    if(match_bits.backend().size() > result.size())
        throw exception("TriggerDescriptorCollection: inconsistent definition of containers");
    size_t n = 0;
    for(; n < match_bits.backend().size(); ++n)
        result[n] = *(match_bits.backend().limbs() + n);
    for(; n < result.size(); ++n)
        result[n] = 0;
    return result;
}

TriggerDescriptorCollection::BitsContainer TriggerDescriptorCollection::ConvertFromRootRepresentation(
    const RootBitsContainer& match_bits)
{
    TriggerDescriptorCollection::BitsContainer result;
    result.backend().resize(match_bits.size(), match_bits.size());
    for(size_t n = 0; n < match_bits.size(); ++n)
        *(result.backend().limbs() + n) = match_bits[n];
    return result;
}

std::shared_ptr<TriggerDescriptorCollection> TriggerDescriptorCollection::Load(const std::string& cfg_name, const Channel& channel)
{
    std::shared_ptr<TriggerDescriptorCollection> triggerDescriptors = std::make_shared<TriggerDescriptorCollection>();
    trigger_tools::SetupDescriptor setup;

    trigger_tools::TriggerFileDescriptorCollection trigger_file_descriptors;
    analysis::ConfigReader config_reader;
    trigger_tools::TriggerFileConfigEntryReader trigger_entry_reader(trigger_file_descriptors);
    config_reader.AddEntryReader("PATTERN", trigger_entry_reader, true);

    trigger_tools::SetupDescriptorCollection setup_file_descriptors;
    trigger_tools::SetupConfigEntryReader setup_entry_reader(setup_file_descriptors);
    config_reader.AddEntryReader("SETUP", setup_entry_reader, false);

    //const std::string triggerCfg_full = edm::FileInPath(cfg_name).fullPath();
    config_reader.ReadConfig(cfg_name);

    if(setup_file_descriptors.size() != 1)
        throw exception("More than 1 setup in Reading Trigger Tools cfg");
    setup = setup_file_descriptors.begin()->second;

    for(const auto& entry : trigger_file_descriptors.get_ordered_by_insertion()) {
        const trigger_tools::TriggerFileDescriptor* trigger_file_descriptor = entry.second;
        if(!trigger_file_descriptor->channels.count(channel)) continue;
        const auto& legs = trigger_file_descriptor->legs;
        std::vector<TriggerDescriptorCollection::Leg> legs_vector;
        for (size_t n = 0; n < legs.size(); ++n){
            const analysis::PropertyList leg_list = analysis::Parse<analysis::PropertyList>(legs.at(n));
            const analysis::LegType type = leg_list.Get<analysis::LegType>("type");
            const double pt = leg_list.Get<double>("pt");
            const double delta_pt = setup.deltaPt_map.at(type);
            boost::optional<double> eta;
            if(leg_list.Has("eta"))
                eta = leg_list.Get<double>("eta");
            bool applyL1match = false;
            if(leg_list.Has("applyL1match"))
                applyL1match = leg_list.Get<bool>("applyL1match");
            const TriggerDescriptorCollection::FilterVector filters = leg_list.GetList<std::string>("filters", false);
            legs_vector.emplace_back(type,pt,delta_pt,eta,applyL1match,filters);
        }
        (*triggerDescriptors).Add(entry.first, legs_vector);
    }

    return triggerDescriptors;
}

const std::vector<std::string>& TriggerDescriptorCollection::GetJetFilters() const { return jet_filters; }

TriggerResults::BitsContainer TriggerResults::GetAcceptBits() const { return accept_bits.to_ullong(); }
TriggerResults::BitsContainer TriggerResults::GetMatchBits() const { return match_bits.to_ullong(); }

void TriggerResults::SetAcceptBits(BitsContainer _accept_bits) { accept_bits = Bits(_accept_bits); }
void TriggerResults::SetMatchBits(BitsContainer _match_bits) { match_bits = Bits(_match_bits); }
void TriggerResults::SetDescriptors(DescriptorsPtr _triggerDescriptors) { triggerDescriptors = _triggerDescriptors; }

bool TriggerResults::Accept(size_t index) const { CheckIndex(index); return accept_bits[index]; }
bool TriggerResults::Match(size_t index) const { CheckIndex(index); return match_bits[index]; }
bool TriggerResults::AcceptAndMatch(size_t index) const
{
    CheckIndex(index);
    return accept_bits[index] && match_bits[index];
}

void TriggerResults::SetAccept(size_t index, bool value) { CheckIndex(index); accept_bits[index] = value; }
void TriggerResults::SetMatch(size_t index, bool value) { CheckIndex(index); match_bits[index] = value; }

bool TriggerResults::Accept(const Pattern& pattern) const { return Accept(GetIndex(pattern)); }
bool TriggerResults::Match(const Pattern& pattern) const { return Match(GetIndex(pattern)); }
bool TriggerResults::AcceptAndMatch(const Pattern& pattern) const { return AcceptAndMatch(GetIndex(pattern)); }

bool TriggerResults::AnyAccept() const { return accept_bits.any(); }
bool TriggerResults::AnyMatch() const { return match_bits.any(); }
bool TriggerResults::AnyAcceptAndMatch() const { return (accept_bits & match_bits).any(); }

bool TriggerResults::MatchEx(size_t index, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    if(!Match(index)) return false;

    const auto& desc = GetTriggerDescriptors().at(index);


    const TriggerDescriptorCollection::Leg& first_leg = desc.lepton_legs.at(0);
    if(pt_firstLeg <= first_leg.pt + first_leg.delta_pt) return false;
    if(desc.lepton_legs.size() > 1){
      const TriggerDescriptorCollection::Leg& second_leg = desc.lepton_legs.at(1);
      if(pt_secondLeg <= second_leg.pt + second_leg.delta_pt) return false;
    }
    
    const size_t n_legs = desc.jet_legs.size();
    if(reco_jet_matches.size() < n_legs) return false;

    if(n_legs == 0) return true;

    if(n_legs == 1) {
        bool match_found = false;
        for(size_t reco_jet_index = 0; !match_found && reco_jet_index < reco_jet_matches.size(); ++reco_jet_index) {
            FilterBitsContainer match_bits = std::numeric_limits<FilterBitsContainer>::max();
            for(unsigned filter_index : desc.jet_legs.at(0).jet_filter_indices) {
                match_bits &= TriggerDescriptorCollection::GetJetFilterMatchBits(reco_jet_matches.at(reco_jet_index), filter_index);
            }
            match_found = match_bits != 0;
        }
        return match_found;
    }

    if(n_legs == 2 && reco_jet_matches.size() == 2) {
        bool match_found = false;
        for(size_t flip = 0; !match_found && flip < n_legs; ++flip) {
            std::vector<FilterBitsContainer> match_bits(n_legs, std::numeric_limits<FilterBitsContainer>::max());
            for(size_t n = 0; n < n_legs; ++n) {
                const size_t reco_jet_index = (n + flip) % 2;
                for(unsigned filter_index : desc.jet_legs.at(n).jet_filter_indices) {
                    match_bits.at(n) &= TriggerDescriptorCollection::GetJetFilterMatchBits(reco_jet_matches.at(reco_jet_index), filter_index);
                }
            }
            const std::bitset<TriggerDescriptorCollection::MaxNumberOfJetFilters> cmb(match_bits.at(0) | match_bits.at(1));
            match_found = match_bits.at(0) != 0 && match_bits.at(1) != 0 && cmb.count() >= n_legs;
        }
        return match_found;
    }

    throw exception("Unsupported number of jet trigger legs.");
}

bool TriggerResults::AcceptAndMatchEx(size_t index, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    if(!AcceptAndMatch(index)) return false;
    return MatchEx(index,pt_firstLeg,pt_secondLeg,reco_jet_matches);
}

bool TriggerResults::MatchEx(const Pattern& pattern, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    return MatchEx(GetIndex(pattern),pt_firstLeg,pt_secondLeg, reco_jet_matches);
}

bool TriggerResults::AcceptAndMatchEx(const Pattern& pattern, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    return AcceptAndMatchEx(GetIndex(pattern), pt_firstLeg, pt_secondLeg, reco_jet_matches);
}

bool TriggerResults::AnyMatchEx(double pt_firstLeg, double pt_secondLeg,const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    if(!AnyMatch()) return false;
    for(size_t n = 0; n < GetTriggerDescriptors().size(); ++n) {
        if(MatchEx(n,pt_firstLeg,pt_secondLeg,reco_jet_matches)) return true;
    }
    return false;
}

bool TriggerResults::AnyAcceptAndMatchEx(double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches) const
{
    if(!AnyAcceptAndMatch()) return false;
    for(size_t n = 0; n < GetTriggerDescriptors().size(); ++n) {
        if(AcceptAndMatchEx(n,pt_firstLeg,pt_secondLeg,reco_jet_matches)) return true;
    }

    return false;
}

void TriggerResults::CheckIndex(size_t index) const
{
    if(index >= MaxNumberOfTriggers)
        throw exception("Trigger index is out of range.");
}

const TriggerDescriptorCollection& TriggerResults::GetTriggerDescriptors() const
{
    if(!triggerDescriptors)
        throw exception("Trigger descriptors not set.");
    return *triggerDescriptors;
}

size_t TriggerResults::GetIndex(const Pattern& pattern) const { return GetTriggerDescriptors().GetIndex(pattern); }

} // namespace nutple
