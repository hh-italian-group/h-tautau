/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "EventTuple.h"

#define SUMMARY_DATA() \
    /* Run statistics */ \
    VAR(UInt_t, exeTime) \
    VAR(ULong64_t, numberOfProcessedEvents) \
    VAR(Double_t, totalShapeWeight) \
    VAR(Double_t, totalShapeWeight_withTopPt) \
    /* Tau ID information */ \
    VAR(std::vector<std::string>, tauId_names) \
    VAR(std::vector<uint32_t>, tauId_keys) \
    /* Trigger information */ \
    VAR(std::vector<Int_t>, triggers_channel) \
    VAR(std::vector<UInt_t>, triggers_index) \
    VAR(std::vector<std::string>, triggers_pattern) \
    VAR(std::vector<UInt_t>, triggers_n_legs) \
    VAR(std::vector<Int_t>, triggerFilters_channel) \
    VAR(std::vector<UInt_t>, triggerFilters_triggerIndex) \
    VAR(std::vector<UInt_t>, triggerFilters_LegId) \
    VAR(std::vector<std::string>, triggerFilters_name) \
    /* MC truth event splitting */ \
    VAR(std::vector<UInt_t>, lhe_n_partons) \
    VAR(std::vector<UInt_t>, lhe_n_b_partons) \
    VAR(std::vector<UInt_t>, lhe_ht10_bin) \
    VAR(std::vector<ULong64_t>, lhe_n_events) \
    /* Skimmer Variables */\
    VAR(std::vector<UInt_t>, file_desc_id) /* vector of File id in TupleSkimmer. */ \
    VAR(std::vector<std::string>, file_desc_name) /* vector of File name in TupleSkimmer. */ \
    VAR(UInt_t, n_splits) /* Number of splits for a file in TupleSkimmer. */ \
    VAR(UInt_t, split_seed) /* Seed for splitting in TupleSkimmer. */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ProdSummary, SummaryTuple, SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, SummaryTuple, SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA


#define EVENT_EXPRESS_DATA() \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, gen_top_pt) /* pt of gen ME top */ \
    VAR(Float_t, gen_topBar_pt) /* pt of gen ME anti-top */ \
    VAR(Float_t, lhe_H_m) /* mass of lhe H */ \
    VAR(Float_t, lhe_hh_m) /* mass of lhe hh pair */ \
    VAR(Float_t, lhe_hh_cosTheta) /* cos(theta) between h and z-axis in the hh reference frame */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, ExpressEvent, ExpressTuple, EVENT_EXPRESS_DATA, "all_events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, ExpressTuple, EVENT_EXPRESS_DATA)
#undef VAR
#undef SUMMARY_DATA

namespace ntuple {
struct GenId {
    size_t n_partons;
    size_t n_b_partons;
    size_t ht10_bin;

    GenId() : n_partons(0), n_b_partons(0), ht10_bin(0) {}
    GenId(size_t _n_partons, size_t _n_b_partons, size_t _ht10_bin) :
        n_partons(_n_partons), n_b_partons(_n_b_partons), ht10_bin(_ht10_bin) {}

    bool operator<(const GenId& other) const
    {
        if(n_partons != other.n_partons) return n_partons < other.n_partons;
        if(n_b_partons != other.n_b_partons) return n_b_partons < other.n_b_partons;
        return ht10_bin < other.ht10_bin;
    }
};

using GenEventCountMap = std::map<GenId, size_t>;

inline GenEventCountMap ExtractGenEventCountMap(const ProdSummary& s)
{
    GenEventCountMap m;
    for(size_t n = 0; n < s.lhe_n_partons.size(); ++n) {
        const GenId id(s.lhe_n_partons.at(n), s.lhe_n_b_partons.at(n), s.lhe_ht10_bin.at(n));
        if(m.count(id))
            throw analysis::exception("Duplicated LHE gen id in prod summary.");
        m[id] = s.lhe_n_events.at(n);
    }
    return m;
}

inline void ConvertGenEventCountMap(ProdSummary& s, const GenEventCountMap& genCountMap)
{
    s.lhe_n_partons.clear();
    s.lhe_n_b_partons.clear();
    s.lhe_ht10_bin.clear();
    s.lhe_n_events.clear();

    for(const auto& bin : genCountMap) {
        s.lhe_n_partons.push_back(static_cast<UInt_t>(bin.first.n_partons));
        s.lhe_n_b_partons.push_back(static_cast<UInt_t>(bin.first.n_b_partons));
        s.lhe_ht10_bin.push_back(static_cast<UInt_t>(bin.first.ht10_bin));
        s.lhe_n_events.push_back(bin.second);
    }
}

inline std::shared_ptr<SummaryTuple> CreateSummaryTuple(const std::string& name, TDirectory* directory,
                                                        bool readMode, TreeState treeState,
                                                        bool ignore_trigger_branches = false)
{
    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        { TreeState::Full, { "totalShapeWeight", "totalShapeWeight_withTopPt",
                             "file_desc_id", "file_desc_name", "n_splits", "split_seed" } },
        { TreeState::Skimmed, { } }
    };

    static const std::set<std::string> trigger_branches = {
        "triggers_channel", "triggers_index", "triggers_pattern", "triggers_n_legs", "triggerFilters_channel",
        "triggerFilters_triggerIndex", "triggerFilters_LegId", "triggerFilters_name"
    };
    auto disabled = disabled_branches.at(treeState);
    if(ignore_trigger_branches)
        disabled.insert(trigger_branches.begin(), trigger_branches.end());

    return std::make_shared<SummaryTuple>(name, directory, readMode, disabled);
}

inline void CheckProdSummaryConsistency(const ProdSummary& s)
{
    if(s.tauId_keys.size() != s.tauId_names.size())
        throw analysis::exception("Inconsistent tauId info in prod summary.");
    const size_t n_trig = s.triggers_channel.size();
    if(s.triggers_index.size() != n_trig || s.triggers_n_legs.size() != n_trig || s.triggers_pattern.size() != n_trig)
        throw analysis::exception("Inconsistent trigger info in prod summary.");
    const size_t n_filter = s.triggerFilters_channel.size();
    if(s.triggerFilters_LegId.size() != n_filter || s.triggerFilters_name.size() != n_filter
            || s.triggerFilters_triggerIndex.size() != n_filter)
        throw analysis::exception("Inconsistent trigger filters info in prod summary.");
    const size_t n_lhe = s.lhe_n_partons.size();
    if(s.lhe_ht10_bin.size() != n_lhe || s.lhe_n_b_partons.size() != n_lhe || s.lhe_n_events.size() != n_lhe)
        throw analysis::exception("Inconsistent LHE info in prod summary.");
}

inline bool CheckProdSummaryCompatibility(const ProdSummary& s1, const ProdSummary& s2, std::ostream* os = nullptr)
{
    if(s1.tauId_keys.size() != s2.tauId_keys.size()) {
        if(os) *os << "Number of tou id keys are not compatible: " << s1.tauId_keys.size() << "!="
                   << s2.tauId_keys.size() << "." << std::endl;
        return false;
    }
    for(size_t n = 0; n < s1.tauId_keys.size(); ++n) {
        if(s1.tauId_keys.at(n) != s2.tauId_keys.at(n) || s1.tauId_names.at(n) != s2.tauId_names.at(n)) {
            if(os) *os << "Tau id key is not compatible: (" << s1.tauId_keys.at(n) << ", "
                       << s1.tauId_names.at(n) << ") != (" << s2.tauId_keys.at(n) << ", "
                       << s1.tauId_names.at(n) << ")." << std::endl;
            return false;
        }
    }
    if(s1.triggers_channel.size() != s2.triggers_channel.size()
            || s1.triggerFilters_channel.size() != s2.triggerFilters_channel.size())
        return false;
    for(size_t n = 0; n < s1.triggers_channel.size(); ++n) {
        if(s1.triggers_channel.at(n) != s2.triggers_channel.at(n)
                || s1.triggers_index.at(n) != s2.triggers_index.at(n)
                || s1.triggers_n_legs.at(n) != s2.triggers_n_legs.at(n)
                || s1.triggers_pattern.at(n) != s2.triggers_pattern.at(n))
            return false;
    }
    for(size_t n = 0; n < s1.triggerFilters_channel.size(); ++n) {
        if(s1.triggerFilters_channel.at(n) != s2.triggerFilters_channel.at(n)
                || s1.triggerFilters_LegId.at(n) != s2.triggerFilters_LegId.at(n)
                || s1.triggerFilters_name.at(n) != s2.triggerFilters_name.at(n)
                || s1.triggerFilters_triggerIndex.at(n) != s2.triggerFilters_triggerIndex.at(n))
            return false;
    }
    return true;
}

inline void MergeProdSummaries(ProdSummary& summary, const ProdSummary& otherSummary)
{
    CheckProdSummaryConsistency(summary);
    CheckProdSummaryConsistency(otherSummary);
    auto genCountMap = ExtractGenEventCountMap(summary);
    if(!CheckProdSummaryCompatibility(summary, otherSummary, &std::cerr))
        throw analysis::exception("Can't merge two incompatible prod summaries.");

    summary.exeTime += otherSummary.exeTime;
    summary.numberOfProcessedEvents += otherSummary.numberOfProcessedEvents;
    summary.totalShapeWeight += otherSummary.totalShapeWeight;
    summary.totalShapeWeight_withTopPt += otherSummary.totalShapeWeight_withTopPt;
    auto otherGenCountMap = ExtractGenEventCountMap(otherSummary);
    for(const auto& bin : otherGenCountMap)
        genCountMap[bin.first] += bin.second;

    ConvertGenEventCountMap(summary, genCountMap);
}

inline ProdSummary MergeSummaryTuple(SummaryTuple& tuple)
{
    ProdSummary summary;
    const Long64_t n_entries = tuple.GetEntries();
    if(n_entries < 1)
        throw analysis::exception("Summary tuple is empty.");
    tuple.GetEntry(0);
    summary = tuple.data();
    CheckProdSummaryConsistency(summary);
    auto genCountMap = ExtractGenEventCountMap(summary);
    for(Long64_t n = 1; n < n_entries; ++n) {
        tuple.GetEntry(n);
        const ProdSummary entry = tuple.data();
        CheckProdSummaryConsistency(entry);
//        try {
            if(!CheckProdSummaryCompatibility(summary, entry, &std::cerr))
                throw analysis::exception("Incompatible prod summaries inside one summary tuple."
                                          " Entry id = %1% out of %2%.") % n % n_entries;
//        } catch(analysis::exception& e) {
//            std::cerr << "ERROR: " << e.what() << std::endl;
//        }

        summary.exeTime += entry.exeTime;
        summary.numberOfProcessedEvents += entry.numberOfProcessedEvents;
        summary.totalShapeWeight += entry.totalShapeWeight;
        summary.totalShapeWeight_withTopPt += entry.totalShapeWeight_withTopPt;
        auto otherGenCountMap = ExtractGenEventCountMap(entry);
        for(const auto& bin : otherGenCountMap)
            genCountMap[bin.first] += bin.second;
    }

    ConvertGenEventCountMap(summary, genCountMap);
    return summary;
}

} // namespace ntuple
