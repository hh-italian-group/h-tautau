/*! Definition of class to represent trigger results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include <vector>
#include <cstdint>
#include <array>
#include <boost/regex.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTypes.h"
#include "h-tautau/Core/include/TriggerFileDescriptor.h"
#include "h-tautau/Core/include/TriggerFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"

namespace analysis {

class TriggerDescriptorCollection {
public:
    using Pattern = std::string;
    using Filter = std::string;
    using BitsContainer = boost::multiprecision::uint256_t;
    using FilterBitsContainer = uint16_t;
    static constexpr size_t MaxNumberOfJetFilters = std::numeric_limits<FilterBitsContainer>::digits;
    static constexpr size_t MaxNumberOfTriggerJets = std::numeric_limits<BitsContainer>::digits
                                                   / MaxNumberOfJetFilters;
    using RootBitsContainerUnit = uint64_t;
    static constexpr size_t NumberOfRootBitsContainerUnits = std::numeric_limits<BitsContainer>::digits
                                                           / std::numeric_limits<RootBitsContainerUnit>::digits;
    using RootBitsContainer = std::array<RootBitsContainerUnit, NumberOfRootBitsContainerUnits>;

    static_assert(std::numeric_limits<RootBitsContainerUnit>::digits
                  == std::numeric_limits<boost::multiprecision::limb_type>::digits,
                  "TriggerDescriptorCollection: inconsistent definition of containers");

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
        double delta_pt;
        boost::optional<double> eta;
        boost::optional<unsigned> run_switch;
        bool applyL1match;
        FilterVector filters;
        boost::optional<FilterVector> legacy_filters;
        std::vector<unsigned> jet_filter_indices;
        Leg(const LegType _type, double _pt, double _delta_pt, boost::optional<double> _eta,
            boost::optional<unsigned> _run_switch, bool _applyL1match,
            const FilterVector& _filters, boost::optional<FilterVector> _legacy_filters);
    };

    struct TriggerDescriptor {
        Pattern pattern;
        std::vector<Leg> lepton_legs, jet_legs;
        boost::regex regex;

        TriggerDescriptor(const Pattern _pattern, const std::vector<Leg>& legs_info);
        bool PatternMatch(const std::string& path_name) const;
        bool RequiresJetMatching() const;
    };

    static FilterBitsContainer GetJetFilterMatchBits(BitsContainer match_bits, unsigned filter_index);
    static RootBitsContainer ConvertToRootRepresentation(BitsContainer match_bits);
    static BitsContainer ConvertFromRootRepresentation(const RootBitsContainer& match_bits);
    static std::shared_ptr<TriggerDescriptorCollection> Load(const std::string& cfg_name, const Channel& channel);

    size_t size() const;
    const TriggerDescriptor& at(size_t index) const;
    const TriggerDescriptor& at(const Pattern& pattern) const;
    void Add(const Pattern& pattern, const std::vector<Leg>& legs);
    bool FindPatternMatch(const std::string& path_name, size_t& index) const;
    size_t GetIndex(const Pattern& pattern) const;

    const std::vector<std::string>& GetJetFilters() const;

private:
    std::vector<TriggerDescriptor> descriptors;
    std::unordered_map<Pattern, size_t> desc_indices;
    mutable std::unordered_map<std::string, size_t> path_index_cache;
    std::vector<std::string> jet_filters;
};

class TriggerResults {
public:
    using BitsContainer = unsigned long long;
    using JetBitsContainer = TriggerDescriptorCollection::BitsContainer;
    using FilterBitsContainer = TriggerDescriptorCollection::FilterBitsContainer;
    static constexpr size_t MaxNumberOfTriggers = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfTriggers>;
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

    template<typename PatternCollection>
    bool AnyAccept(const PatternCollection& patterns) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return Accept(pattern); });
    }

    bool AnyAccept() const;
    bool AnyMatch() const;
    bool AnyAcceptAndMatch() const;

    bool MatchEx(size_t index, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;
    bool AcceptAndMatchEx(size_t index, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;
    bool MatchEx(const Pattern& pattern, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;
    bool AcceptAndMatchEx(const Pattern& pattern, double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;

    template<typename PatternCollection>
    bool AnyAcceptAndMatchEx(const PatternCollection& patterns, double pt_firstLeg, double pt_secondLeg,
                             const std::vector<JetBitsContainer>& reco_jet_matches = {}) const
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&](const Pattern& pattern) { return AcceptAndMatchEx(pattern, pt_firstLeg, pt_secondLeg, reco_jet_matches); });
    }

    bool AnyMatchEx(double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;
    bool AnyAcceptAndMatchEx(double pt_firstLeg, double pt_secondLeg, const std::vector<JetBitsContainer>& reco_jet_matches = {}) const;

private:
    void CheckIndex(size_t index) const;
    const TriggerDescriptorCollection& GetTriggerDescriptors() const;
    size_t GetIndex(const Pattern& pattern) const;

private:
    Bits accept_bits, match_bits;
    DescriptorsPtr triggerDescriptors;
};

} // namespace nutple
