/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/EventTuple.h"
#include "AnalysisTools/Core/include/Tools.h"

namespace ntuple {

const LegPair LegPair::Undefined(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());

LegPair::LegPair() : std::pair<size_t, size_t>(Undefined) {}
LegPair::LegPair(size_t _first, size_t _second) : std::pair<size_t, size_t>(_first, _second) {}
LegPair::LegPair(const std::pair<size_t, size_t>& p) : std::pair<size_t, size_t>(p) {}

size_t LegPair::Get(size_t position) const
{
    if(position == 1) return first;
    if(position == 2) return second;
    throw analysis::exception("LegPair: position is out of range");
}

size_t LegPair::ToIndex() const { return first * 1000 + second; }
bool LegPair::IsDefined() const { return first != Undefined.first && second != Undefined.second; }
bool LegPair::Contains(size_t i) const { return first == i || second == i; }

LegPair LegPair::FromIndex(size_t index)
{
    LegPair pair;
    pair.second = index % 1000;
    pair.first = (index - pair.second) / 1000;
    return pair;
}

std::shared_ptr<EventTuple> CreateEventTuple(const std::string& name, TDirectory* directory,
                                                    bool readMode, TreeState treeState)
{
    static const std::set<std::string> weight_branches = {
        "weight_pu", "weight_pu_up", "weight_pu_down", "weight_dy", "weight_ttbar", "weight_wjets",
        "weight_bsm_to_sm", "weight_top_pt", "weight_xs", "weight_xs_withTopPt", "weight_total", "weight_total_withTopPt"
    };

    static const std::set<std::string> gen_study_branches = {
        "sample_type", "sample_year", "mass_point", "spin", "node",
    };

    static const std::set<std::string> skimmer_branches = {
        "period", "isData", "file_desc_id", "split_id",
    };

    static const std::set<std::string> SVfit_branches = {
        "SVfit_htt_index", "SVfit_is_valid", "SVfit_p4", "SVfit_p4_error", "SVfit_mt", "SVfit_mt_error",
        "SVfit_unc_source", "SVfit_unc_scale",
    };

    static const std::set<std::string> kinFit_branches = {
        "kinFit_htt_index", "kinFit_hbb_index", "kinFit_m", "kinFit_chi2", "kinFit_convergence",
        "kinFit_unc_source", "kinFit_unc_scale",
    };

    static const std::set<std::string> HHbtag_branches = {
        "jet_HHbtag_htt_index", "jet_HHbtag_jet_index", "jet_HHbtag_unc_source", "jet_HHbtag_unc_scale",
        "jet_HHbtag_value",
    };

    static const std::map<TreeState, std::set<std::string>> disabled_branches = {
        {
            TreeState::Full,
            analysis::tools::union_sets({
                weight_branches, gen_study_branches, skimmer_branches, SVfit_branches, kinFit_branches,
                HHbtag_branches,
            })
        },
        { TreeState::Skimmed, { } },
    };

    const auto& disabled = disabled_branches.at(treeState);
    return std::make_shared<EventTuple>(name, directory, readMode, disabled);
}

} // namespace ntuple
