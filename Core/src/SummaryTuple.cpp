/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/SummaryTuple.h"

namespace ntuple {
GenId::GenId() : n_partons(0), n_b_partons(0), ht10_bin(0) {}
GenId::GenId(size_t _n_partons, size_t _n_b_partons, size_t _ht10_bin) :
    n_partons(_n_partons), n_b_partons(_n_b_partons), ht10_bin(_ht10_bin) {}

bool GenId::operator<(const GenId& other) const
{
    if(n_partons != other.n_partons) return n_partons < other.n_partons;
    if(n_b_partons != other.n_b_partons) return n_b_partons < other.n_b_partons;
    return ht10_bin < other.ht10_bin;
}

GenEventCountMap ExtractGenEventCountMap(const ProdSummary& s)
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

void ConvertGenEventCountMap(ProdSummary& s, const GenEventCountMap& genCountMap)
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

GenEventTypeCountMap ExtractGenEventTypeCountMap(const ProdSummary& s)
{
    GenEventTypeCountMap m;
    for(size_t n = 0; n < s.genEventType.size(); ++n) {
        analysis::GenEventType genEventType = static_cast<analysis::GenEventType>(s.genEventType.at(n));
        if(m.count(genEventType))
            throw analysis::exception("Duplicated genEventType in prod summary.");
        m[genEventType] = s.genEventType_n_events.at(n);
    }
    return m;
}

void ConvertGenEventTypeCountMap(ProdSummary& s, const GenEventTypeCountMap& genCountMap)
{
    s.genEventType.clear();
    s.genEventType_n_events.clear();

    for(const auto& bin : genCountMap) {
        s.genEventType.push_back(static_cast<Int_t>(bin.first));
        s.genEventType_n_events.push_back(static_cast<ULong64_t>(bin.second));
    }
}


std::shared_ptr<SummaryTuple> CreateSummaryTuple(const std::string& name, TDirectory* directory,
                                                        bool readMode, TreeState treeState)
{
    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        { TreeState::Full, { "totalShapeWeight", "totalShapeWeight_withTopPt",
                             "file_desc_id", "file_desc_name", "n_splits", "split_seed" } },
        { TreeState::Skimmed, { } }
    };

    const auto& disabled = disabled_branches.at(treeState);
    return std::make_shared<SummaryTuple>(name, directory, readMode, disabled);
}

void CheckProdSummaryConsistency(const ProdSummary& s)
{
    if(s.triggers_pattern.size() != s.triggers_channel.size())
        throw analysis::exception("Inconsistent trigger info in prod summary.");
    const size_t n_lhe = s.lhe_n_partons.size();
    if(s.lhe_ht10_bin.size() != n_lhe || s.lhe_n_b_partons.size() != n_lhe || s.lhe_n_events.size() != n_lhe)
        throw analysis::exception("Inconsistent LHE info in prod summary.");
    const size_t n_genEventType = s.genEventType.size();
    if(s.genEventType_n_events.size() != n_genEventType)
        throw analysis::exception("Inconsistent genEventType info in prod summary.");
}

bool CheckProdSummaryCompatibility(const ProdSummary& s1, const ProdSummary& s2)
{
    if(s1.triggers_channel.size() != s2.triggers_channel.size())
        return false;
    for(size_t n = 0; n < s1.triggers_channel.size(); ++n) {
        if(s1.triggers_channel.at(n) != s2.triggers_channel.at(n)
                || s1.triggers_pattern.at(n) != s2.triggers_pattern.at(n))
            return false;
    }
    return true;
}

void MergeProdSummaries(ProdSummary& summary, const ProdSummary& otherSummary)
{
    CheckProdSummaryConsistency(summary);
    CheckProdSummaryConsistency(otherSummary);
    auto genCountMap = ExtractGenEventCountMap(summary);
    auto genEventTypeCountMap = ExtractGenEventTypeCountMap(summary);
    if(!CheckProdSummaryCompatibility(summary, otherSummary))
        throw analysis::exception("Can't merge two incompatible prod summaries.");

    summary.exeTime += otherSummary.exeTime;
    summary.numberOfProcessedEvents += otherSummary.numberOfProcessedEvents;
    summary.totalShapeWeight += otherSummary.totalShapeWeight;
    summary.totalShapeWeight_withTopPt += otherSummary.totalShapeWeight_withTopPt;
    auto otherGenCountMap = ExtractGenEventCountMap(otherSummary);
    auto otherGenEventTypeCountMap = ExtractGenEventTypeCountMap(otherSummary);
    for(const auto& bin : otherGenCountMap)
        genCountMap[bin.first] += bin.second;

    for(const auto& bin : otherGenEventTypeCountMap)
        genEventTypeCountMap[bin.first] += bin.second;

    ConvertGenEventCountMap(summary, genCountMap);
    ConvertGenEventTypeCountMap(summary, genEventTypeCountMap);
}

ProdSummary MergeSummaryTuple(SummaryTuple& tuple)
{
    ProdSummary summary;
    const Long64_t n_entries = tuple.GetEntries();
    if(n_entries < 1)
        throw analysis::exception("Summary tuple is empty.");
    tuple.GetEntry(0);
    summary = tuple.data();
    CheckProdSummaryConsistency(summary);
    auto genCountMap = ExtractGenEventCountMap(summary);
    auto genEventTypeCountMap = ExtractGenEventTypeCountMap(summary);
    for(Long64_t n = 1; n < n_entries; ++n) {
        tuple.GetEntry(n);
        const ProdSummary entry = tuple.data();
        CheckProdSummaryConsistency(entry);
//        try {
            if(!CheckProdSummaryCompatibility(summary, entry))
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
        auto otherGenEventTypeCountMap = ExtractGenEventTypeCountMap(entry);
        for(const auto& bin : otherGenCountMap)
            genCountMap[bin.first] += bin.second;

        for(const auto& bin : otherGenEventTypeCountMap)
            genEventTypeCountMap[bin.first] += bin.second;
    }

    ConvertGenEventCountMap(summary, genCountMap);
    ConvertGenEventTypeCountMap(summary, genEventTypeCountMap);
    return summary;
}


std::shared_ptr<ExpressTuple> CreateExpressTuple(const std::string& name, TDirectory* directory,
                                                        bool readMode, TreeState treeState)
{
    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        { TreeState::Full, { "file_desc_id", "split_id"} },
        { TreeState::Skimmed, {} }
    };

    const auto& disabled = disabled_branches.at(treeState);
    return std::make_shared<ExpressTuple>(name, directory, readMode, disabled);
}

} // namespace ntuple
