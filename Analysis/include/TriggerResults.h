/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include <vector>
#include <boost/regex.hpp>
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTypes.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

class TriggerDescriptorCollection {
public:
    using Pattern = std::string;
    using Filter = std::string;
    using BitsContainer = unsigned long long;
    static constexpr size_t MaxNumberOfTriggers = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfTriggers>;
    static constexpr size_t MaxNumberOfJetFilters = 4;
    static constexpr size_t MaxNumberOfTriggerJets = std::numeric_limits<BitsContainer>::digits
                                                   / MaxNumberOfJetFilters;

    using PatternContainer = std::vector<Pattern>;
    using FilterVector = std::vector<Filter>;
    using FilterContainer = std::map<size_t, FilterVector>;

    struct JetTriggerObjectCollection {
        BitsContainer match_bits;
        std::vector<LorentzVector> momentums;

        JetTriggerObjectCollection() : match_bits(0) {}

        template<typename LVector>
        size_t Add(const LVector& jet_momentum)
        {
            if(momentums.size() == MaxNumberOfTriggerJets)
                throw exception("Maximal number of trigger jets is reached.");
            momentums.emplace_back(jet_momentum);
            return momentums.size() - 1;
        }

        bool GetJetFilterMatchBit(size_t filter_index, size_t jet_index) const
        {
            return ( match_bits >> (filter_index * MaxNumberOfTriggerJets + jet_index) ) & 1;
        }

        void SetJetFilterMatchBit(size_t filter_index, size_t jet_index, bool match_result)
        {
            const BitsContainer mask = BitsContainer(1) << (filter_index * MaxNumberOfTriggerJets + jet_index);
            match_bits = match_result ? match_bits | mask : match_bits & ~mask;
        }
    };

    struct Leg {
        LegType type;
        double pt;
        FilterVector filters;
        std::vector<unsigned> jet_filter_indices;

        Leg(const LegType _type, double _pt, const FilterVector& _filters)
            : type(_type), pt(_pt), filters(_filters) { }
    };

    struct TriggerDescriptor {
        Pattern pattern;
        std::vector<Leg> lepton_legs, jet_legs;
        boost::regex regex;

        TriggerDescriptor(const Pattern _pattern, const std::vector<Leg>& legs_info)
            : pattern(_pattern)
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

        bool PatternMatch(const std::string& path_name) const { return boost::regex_match(path_name, regex); }
        bool RequiresJetMatching() const { return !jet_legs.empty(); }
    };

    size_t size() const { return descriptors.size(); }

    const TriggerDescriptor& at(size_t index) const
    {
        if(index >= descriptors.size())
            throw exception("Index out of trigger descriptors range.");
        return descriptors.at(index);
    }

    const TriggerDescriptor& at(const Pattern& pattern) const
    {
        size_t index = GetIndex(pattern);
        return descriptors.at(index);
    }

    void Add(const Pattern& pattern, const std::vector<Leg>& legs)
    {
        if(desc_indices.count(pattern))
            throw exception("Duplicated trigger pattern '%1%'.") % pattern;
        if(descriptors.size() == MaxNumberOfTriggers)
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

    bool FindPatternMatch(const std::string& path_name, size_t& index)
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

    size_t GetIndex(const Pattern& pattern) const
    {
        if(!desc_indices.count(pattern))
            throw exception("Trigger '%1%' not found.") % pattern;
        return desc_indices.at(pattern);
    }

    BitsContainer GetJetFilterMatchBits(BitsContainer match_bits, unsigned filter_index) const
    {
        static constexpr BitsContainer mask((1 << MaxNumberOfTriggerJets) - 1);
        return ( match_bits >> (filter_index * MaxNumberOfTriggerJets) ) & mask;
    }

    const std::vector<std::string>& GetJetFilters() const { return jet_filters; }

private:
    std::vector<TriggerDescriptor> descriptors;
    std::unordered_map<Pattern, size_t> desc_indices;
    std::unordered_map<std::string, size_t> path_index_cache;
    std::vector<std::string> jet_filters;
};

class TriggerResults {
public:
    using BitsContainer = TriggerDescriptorCollection::BitsContainer;
    static constexpr size_t MaxNumberOfTriggers = TriggerDescriptorCollection::MaxNumberOfTriggers;
    static constexpr size_t MaxNumberOfTriggerJets = TriggerDescriptorCollection::MaxNumberOfTriggerJets;
    using Bits = TriggerDescriptorCollection::Bits;
    using Pattern = TriggerDescriptorCollection::Pattern;

    using DescriptorsPtr = std::shared_ptr<const TriggerDescriptorCollection>;

    BitsContainer GetAcceptBits() const { return accept_bits.to_ullong(); }
    BitsContainer GetMatchBits() const { return match_bits.to_ullong(); }

    void SetAcceptBits(BitsContainer _accept_bits) { accept_bits = Bits(_accept_bits); }
    void SetMatchBits(BitsContainer _match_bits) { match_bits = Bits(_match_bits); }
    void SetDescriptors(DescriptorsPtr _triggerDescriptors) { triggerDescriptors = _triggerDescriptors; }

    bool Accept(size_t index) const { CheckIndex(index); return accept_bits[index]; }
    bool Match(size_t index) const { CheckIndex(index); return match_bits[index]; }
    bool AcceptAndMatch(size_t index) const { CheckIndex(index); return accept_bits[index] && match_bits[index]; }
    void SetAccept(size_t index, bool value) { CheckIndex(index); accept_bits[index] = value; }
    void SetMatch(size_t index, bool value) { CheckIndex(index); match_bits[index] = value; }

    bool Accept(const Pattern& pattern) const { return Accept(GetIndex(pattern)); }
    bool Match(const Pattern& pattern) const { return Match(GetIndex(pattern)); }
    bool AcceptAndMatch(const Pattern& pattern) const { return AcceptAndMatch(GetIndex(pattern)); }

    template<typename PatternCollection>
    bool AnyAcceptAndMatch(const PatternCollection& patterns) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return AcceptAndMatch(pattern); });
    }

    bool AnyAccpet() const { return accept_bits.any(); }
    bool AnyMatch() const { return match_bits.any(); }
    bool AnyAcceptAndMatch() const { return (accept_bits & match_bits).any(); }

    bool MatchEx(size_t index, const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        if(!Match(index)) return false;

        const auto& desc = GetTriggerDescriptors().at(index);
        const size_t n_legs = desc.jet_legs.size();
        if(reco_jet_matches.size() < n_legs) return false;

        if(n_legs == 0) return true;

        if(n_legs == 1) {
            bool match_found = false;
            for(size_t reco_jet_index = 0; !match_found && reco_jet_index < reco_jet_matches.size(); ++reco_jet_index) {
                BitsContainer match_bits = std::numeric_limits<BitsContainer>::max();
                for(unsigned filter_index : desc.jet_legs.at(0).jet_filter_indices) {
                    match_bits &= GetJetFilterMatchBits(reco_jet_matches.at(reco_jet_index), filter_index);
                }
                match_found = match_bits != 0;
            }
            return match_found;
        }

        if(n_legs == 2 && reco_jet_matches.size() == 2) {
            bool match_found = false;
            for(size_t flip = 0; !match_found && flip < n_legs; ++flip) {
                std::vector<BitsContainer> match_bits(n_legs, std::numeric_limits<BitsContainer>::max());
                for(size_t n = 0; n < n_legs; ++n) {
                    const size_t reco_jet_index = (n + flip) % 2;
                    for(unsigned filter_index : desc.jet_legs.at(n).jet_filter_indices) {
                        match_bits.at(n) &= GetJetFilterMatchBits(reco_jet_matches.at(reco_jet_index), filter_index);
                    }
                }
                const Bits cmb(match_bits.at(0) | match_bits.at(1));
                match_found = match_bits.at(0) != 0 && match_bits.at(1) != 0 && cmb.count() >= n_legs;
            }
            return match_found;
        }

        throw exception("Unsupported number of jet trigger legs.");
    }

    bool AcceptAndMatchEx(size_t index, const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        if(!AcceptAndMatch(index)) return false;
        return MatchEx(index, reco_jet_matches);
    }

    bool MatchEx(const Pattern& pattern, const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        return MatchEx(GetIndex(pattern), reco_jet_matches);
    }

    bool AcceptAndMatchEx(const Pattern& pattern, const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        return AcceptAndMatchEx(GetIndex(pattern), reco_jet_matches);
    }

    template<typename PatternCollection>
    bool AnyAcceptAndMatchEx(const PatternCollection& patterns,
                             const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return AcceptAndMatchEx(pattern, reco_jet_matches); });
    }

    bool AnyMatchEx(const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        if(!AnyMatch()) return false;
        for(size_t n = 0; n < GetTriggerDescriptors().size(); ++n) {
            if(MatchEx(n, reco_jet_matches)) return true;
        }
        return false;
    }

    bool AnyAcceptAndMatchEx(const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        if(!AnyAcceptAndMatch()) return false;
        for(size_t n = 0; n < GetTriggerDescriptors().size(); ++n) {
            if(AcceptAndMatchEx(n, reco_jet_matches)) return true;
        }

        return false;
    }

private:
    void CheckIndex(size_t index) const
    {
        if(index >= MaxNumberOfTriggers)
            throw exception("Trigger index is out of range.");
    }

    const TriggerDescriptorCollection& GetTriggerDescriptors() const
    {
        if(!triggerDescriptors)
            throw exception("Trigger descriptors not set.");
        return *triggerDescriptors;
    }

    size_t GetIndex(const Pattern& pattern) const { return GetTriggerDescriptors().GetIndex(pattern); }

    BitsContainer GetJetFilterMatchBits(BitsContainer match_bits, unsigned filter_index) const
    {
        static constexpr BitsContainer mask((1 << MaxNumberOfTriggerJets) - 1);
        return ( match_bits >> (filter_index * MaxNumberOfTriggerJets) ) & mask;
    }


private:
    Bits accept_bits, match_bits;
    DescriptorsPtr triggerDescriptors;
};

} // namespace nutple
