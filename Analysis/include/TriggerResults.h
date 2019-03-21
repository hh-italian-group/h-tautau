/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include <vector>
#include <boost/regex.hpp>

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Core/include/AnalysisTypes.h"

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

        JetTriggerObjectCollection();

        template<typename LVector>
        size_t Add(const LVector& jet_momentum)
        {
            if(momentums.size() == MaxNumberOfTriggerJets)
                throw exception("Maximal number of trigger jets is reached.");
            momentums.emplace_back(jet_momentum);
            return momentums.size() - 1;
        }

        bool GetJetFilterMatchBit(size_t filter_index, size_t jet_index) const;
        void SetJetFilterMatchBit(size_t filter_index, size_t jet_index, bool match_result);
    };

    struct Leg {
        LegType type;
        double pt;
        boost::optional<double> eta;
        FilterVector filters;
        std::vector<unsigned> jet_filter_indices;
        Leg(const LegType _type, double _pt, boost::optional<double> _eta, const FilterVector& _filters);
    };

    struct TriggerDescriptor {
        Pattern pattern;
        std::vector<Leg> lepton_legs, jet_legs;
        boost::regex regex;

        TriggerDescriptor(const Pattern _pattern, const std::vector<Leg>& legs_info);
        bool PatternMatch(const std::string& path_name) const;
        bool RequiresJetMatching() const;
    };

    size_t size() const;
    const TriggerDescriptor& at(size_t index) const;
    const TriggerDescriptor& at(const Pattern& pattern) const;
    void Add(const Pattern& pattern, const std::vector<Leg>& legs);
    bool FindPatternMatch(const std::string& path_name, size_t& index);
    size_t GetIndex(const Pattern& pattern) const;
    BitsContainer GetJetFilterMatchBits(BitsContainer match_bits, unsigned filter_index) const;
    const std::vector<std::string>& GetJetFilters() const;

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

    BitsContainer GetAcceptBits() const;
    BitsContainer GetMatchBits() const;

    void SetAcceptBits(BitsContainer _accept_bits);
    void SetMatchBits(BitsContainer _match_bits);
    void SetDescriptors(DescriptorsPtr _triggerDescriptors);

    bool Accept(size_t index) const;
    bool Match(size_t index) const;
    bool AcceptAndMatch(size_t index) const;
    void SetAccept(size_t index, bool value);
    void SetMatch(size_t index, bool value);

    bool Accept(const Pattern& pattern) const;
    bool Match(const Pattern& pattern) const;
    bool AcceptAndMatch(const Pattern& pattern) const;

    template<typename PatternCollection>
    bool AnyAcceptAndMatch(const PatternCollection& patterns) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return AcceptAndMatch(pattern); });
    }

    bool AnyAccpet() const;
    bool AnyMatch() const;
    bool AnyAcceptAndMatch() const;

    bool MatchEx(size_t index, const std::vector<BitsContainer>& reco_jet_matches = {}) const;
    bool AcceptAndMatchEx(size_t index, const std::vector<BitsContainer>& reco_jet_matches = {}) const;
    bool MatchEx(const Pattern& pattern, const std::vector<BitsContainer>& reco_jet_matches = {}) const;
    bool AcceptAndMatchEx(const Pattern& pattern, const std::vector<BitsContainer>& reco_jet_matches = {}) const;

    template<typename PatternCollection>
    bool AnyAcceptAndMatchEx(const PatternCollection& patterns,
                             const std::vector<BitsContainer>& reco_jet_matches = {}) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return AcceptAndMatchEx(pattern, reco_jet_matches); });
    }

    bool AnyMatchEx(const std::vector<BitsContainer>& reco_jet_matches = {}) const;
    bool AnyAcceptAndMatchEx(const std::vector<BitsContainer>& reco_jet_matches = {}) const;

private:
    void CheckIndex(size_t index) const;
    const TriggerDescriptorCollection& GetTriggerDescriptors() const;
    size_t GetIndex(const Pattern& pattern) const;
    BitsContainer GetJetFilterMatchBits(BitsContainer match_bits, unsigned filter_index) const;

private:
    Bits accept_bits, match_bits;
    DescriptorsPtr triggerDescriptors;
};

} // namespace nutple
