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



void MergeProdSummaries(ProdSummary& summary, const ProdSummary& otherSummary)
{
    summary.exeTime += otherSummary.exeTime;
    summary.numberOfProcessedEvents += otherSummary.numberOfProcessedEvents;
    summary.totalShapeWeight += otherSummary.totalShapeWeight;
    summary.totalShapeWeight_withTopPt += otherSummary.totalShapeWeight_withTopPt;
}

ProdSummary MergeSummaryTuple(SummaryTuple& tuple)
{
    ProdSummary summary;
    const Long64_t n_entries = tuple.GetEntries();
    if(n_entries < 1)
        throw analysis::exception("Summary tuple is empty.");
    tuple.GetEntry(0);
    summary = tuple.data();
    for(Long64_t n = 1; n < n_entries; ++n) {
        tuple.GetEntry(n);
        const ProdSummary entry = tuple.data();

        summary.exeTime += entry.exeTime;
        summary.numberOfProcessedEvents += entry.numberOfProcessedEvents;
        summary.totalShapeWeight += entry.totalShapeWeight;
        summary.totalShapeWeight_withTopPt += entry.totalShapeWeight_withTopPt;
    }

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
