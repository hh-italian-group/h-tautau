/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include <boost/regex.hpp>
#include "AnalysisTools/Core/include/exception.h"

namespace analysis {

class TriggerDescriptors {
public:
    using Pattern = std::string;
    using Filter = std::string;
    using PatternContainer = std::vector<Pattern>;
    using FilterVector = std::vector<Filter>;
    using FilterContainer = std::map<size_t, FilterVector>;

    const PatternContainer& GetPatterns() const { return patterns; }
    size_t size() const { return patterns.size(); }

    const FilterContainer& GetFilters(size_t index) const
    {
        CheckIndex(index);
        return filters.at(index);
    }

    const FilterVector& GetFilters(size_t index, size_t leg_id) const
    {
        static const FilterVector empty = {};
        const auto& f = GetFilters(index);
        return f.count(leg_id) ? f.at(leg_id) : empty;
    }

    size_t GetNumberOfLegs(size_t index) const
    {
        CheckIndex(index);
        return number_of_legs.at(index);
    }

    void Add(const Pattern& pattern, size_t n_legs, const FilterContainer& leg_filters)
    {
        static const std::string regex_format = "^%1%[0-9]+$";
        if(pattern_indices.count(pattern))
            throw exception("Duplicated trigger pattern '%1%'.") % pattern;
        pattern_indices[pattern] = patterns.size();
        patterns.push_back(pattern);
        number_of_legs.push_back(n_legs);
        filters.push_back(leg_filters);
        const std::string regex_str = boost::str(boost::format(regex_format) % pattern);
        regexes.push_back(boost::regex(regex_str));
    }

    bool PatternMatch(const std::string& path_name, size_t index) const
    {
        CheckIndex(index);
        return boost::regex_match(path_name, regexes.at(index));
    }

    bool FindPatternMatch(const std::string& path_name, size_t& index) const
    {
        for(index = 0; index < patterns.size(); ++index)
            if(PatternMatch(path_name, index)) return true;
        return false;
    }

    size_t GetIndex(const Pattern& pattern) const
    {
        if(!pattern_indices.count(pattern))
            throw exception("Trigger '%1%' not found.") % pattern;
        return pattern_indices.at(pattern);
    }

private:
    void CheckIndex(size_t index) const
    {
        if(index >= patterns.size())
            throw exception("Trigger pattern index is out of range.");
    }

private:
    PatternContainer patterns;
    std::vector<size_t> number_of_legs;
    std::vector<FilterContainer> filters;
    std::vector<boost::regex> regexes;
    std::map<Pattern, size_t> pattern_indices;
};

class TriggerResults {
public:
    using BitsContainer = unsigned long long;
    static constexpr size_t MaxNumberOfTriggers = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfTriggers>;
    using DescriptorsPtr = std::shared_ptr<const TriggerDescriptors>;
    using Pattern = TriggerDescriptors::Pattern;

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
