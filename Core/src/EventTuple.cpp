/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/EventTuple.h"

namespace ntuple {


size_t LegPairToIndex(const LegPair& pair)
{
    return pair.first * 1000 + pair.second;
}

LegPair LegIndexToPair(size_t index)
{
    LegPair pair;
    pair.second = index % 1000;
    pair.first = (index - pair.second) / 1000;
    return pair;
}

LegPair UndefinedJetPair()
{
    static LegPair pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    return pair;
}

std::shared_ptr<EventTuple> CreateEventTuple(const std::string& name, TDirectory* directory,
                                                    bool readMode, TreeState treeState)
{
    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        { TreeState::Full, { "n_jets", "ht_other_jets", "weight_pu", "weight_lepton_trig", "weight_lepton_id_iso",
                             "weight_tau_id", "weight_btag", "weight_btag_up", "weight_btag_down", "weight_dy",
                             "weight_ttbar", "weight_wjets", "weight_bsm_to_sm", "weight_top_pt", "weight_xs",
                             "weight_xs_withTopPt", "weight_total", "weight_total_withTopPt", "isData", "file_desc_id",
                             "split_id" } },
        { TreeState::Skimmed, { } }
    };

    const auto& disabled = disabled_branches.at(treeState);
    return std::make_shared<EventTuple>(name, directory, readMode, disabled);
}

} // namespace ntuple
