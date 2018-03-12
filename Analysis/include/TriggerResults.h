/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include <boost/regex.hpp>
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTypes.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

class TriggerDescriptorCollection {
public:
    using Pattern = std::string;
    using Filter = std::string;
    using PatternContainer = std::vector<Pattern>;
    using FilterVector = std::vector<Filter>;
    using FilterContainer = std::map<size_t, FilterVector>;

    struct Leg{
        LegType type;
        double pt;
        FilterVector filters;

        Leg(const LegType _type, const double _pt, const TriggerDescriptorCollection::FilterVector _filters)
            : type(_type), pt(_pt), filters(_filters) { }

    };

    struct TriggerDescriptor{
        Pattern pattern;
        std::vector<Leg> legs_info;
        boost::regex regex;

        TriggerDescriptor(const Pattern _pattern, const std::vector<Leg> _legs_info)
            : pattern(_pattern), legs_info(_legs_info)
        {
            static const std::string regex_format = "^%1%[0-9]+$";
            const std::string regex_str = boost::str(boost::format(regex_format) % pattern);
            regex = boost::regex(regex_str);
        }

        bool PatternMatch(const std::string& path_name) const
        {
            return boost::regex_match(path_name, regex);
        }


    };

    size_t Size() const
    {
        return pattern_structs.size();
    }

    TriggerDescriptorCollection::TriggerDescriptor GetTriggerDescriptor(size_t index) const
    {
        return pattern_structs.at(index);
    }

    TriggerDescriptorCollection::TriggerDescriptor GetTriggerDescriptor(const Pattern& pattern)
    {
        size_t index = GetIndex(pattern);
        return pattern_structs.at(index);
    }
    
    void Add(const Pattern& pattern, std::vector<Leg> legs)
    {
        if(pattern_indices.count(pattern))
            throw exception("Duplicated trigger pattern '%1%'.") % pattern;
        pattern_indices[pattern] = legs.size();
        TriggerDescriptorCollection::TriggerDescriptor pattern_struct(pattern,legs);
        pattern_structs.push_back(pattern_struct);
    }

    bool FindPatternMatch(const std::string& path_name, size_t& index) const
    {
        for(index = 0; index < pattern_structs.size(); ++index){
            const TriggerDescriptorCollection::TriggerDescriptor pattern_struct = pattern_structs.at(index);
            if(pattern_struct.PatternMatch(path_name)) return true;
        }
        return false;
    }


    size_t GetIndex(const Pattern& pattern) const
    {
        if(!pattern_indices.count(pattern))
            throw exception("Trigger '%1%' not found.") % pattern;
        return pattern_indices.at(pattern);
    }

private:
    std::vector<TriggerDescriptorCollection::TriggerDescriptor> pattern_structs;
    std::map<Pattern, size_t> pattern_indices;
};

class TriggerResults {
public:
    using BitsContainer = unsigned long long;
    static constexpr size_t MaxNumberOfTriggers = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfTriggers>;
    using DescriptorsPtr = std::shared_ptr<const TriggerDescriptorCollection>;
    using Pattern = TriggerDescriptorCollection::Pattern;

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

private:
    void CheckIndex(size_t index) const
    {
        if(index >= MaxNumberOfTriggers)
            throw exception("Trigger index is out of range.");
    }

    size_t GetIndex(const Pattern& pattern) const
    {
        if(!triggerDescriptors)
            throw exception("Trigger descriptors not set.");
        return triggerDescriptors->GetIndex(pattern);
    }

private:
    Bits accept_bits, match_bits;
    DescriptorsPtr triggerDescriptors;
};

} // namespace nutple
