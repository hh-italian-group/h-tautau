/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/EventTuple.h"

namespace ntuple {

size_t NumberOfCombinationPairs(size_t n_jets) { return n_jets * (n_jets - 1); }

size_t CombinationPairToIndex(const JetPair& pair, size_t n_jets)
{
    const size_t min = std::min(pair.first, pair.second);
    const size_t max = std::max(pair.first, pair.second);
    if(n_jets < 2 || min == max || max >= n_jets)
        throw analysis::exception("bad combination pair (%1%, %2%) for n b-jets = %3%.")
            % pair.first % pair.second % n_jets;
    size_t index = pair.first * (n_jets - 1) + pair.second;
    if(pair.first < pair.second)
        --index;
    return index;
}

JetPair CombinationIndexToPair(size_t index, size_t n_jets)
{
    if(n_jets < 2 || index >= NumberOfCombinationPairs(n_jets))
        throw analysis::exception("bad combination index = %1% for n b-jets = %2%.") % index % n_jets;

    JetPair pair;
    pair.second = index % (n_jets - 1);
    pair.first = (index - pair.second) / (n_jets - 1);
    if(pair.first <= pair.second)
        ++pair.second;
    return pair;
}

JetPair UndefinedJetPair()
{
    static JetPair pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    return pair;
}

std::shared_ptr<EventTuple> CreateEventTuple(const std::string& name, TDirectory* directory,
                                                    bool readMode, TreeState treeState)
{
    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        { TreeState::Full, { "n_jets", "ht_other_jets", "weight_pu", "weight_lepton_trig", "weight_lepton_id_iso",
                             "weight_tau_id", "weight_btag", "weight_btag_up", "weight_btag_down", "weight_dy",
                             "weight_ttbar", "weight_wjets", "weight_bsm_to_sm", "weight_top_pt", "weight_xs",
                             "weight_xs_withTopPt", "weight_total", "weight_total_withTopPt", "file_desc_id",
                             "split_id", "isData", "SVfit_Higges_indexes" } },
        { TreeState::Skimmed, { } }
    };

    const auto& disabled = disabled_branches.at(treeState);
    return std::make_shared<EventTuple>(name, directory, readMode, disabled);
}

} // namespace ntuple
