/*! Common simple types for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/AnalysisTypes.h"

namespace analysis {

const std::set<UncertaintyScale>& GetAllUncertaintyScales()
{
    static const std::set<UncertaintyScale> all_scales(
        EnumNameMap<UncertaintyScale>::GetDefault().GetEnumEntries().begin(),
        EnumNameMap<UncertaintyScale>::GetDefault().GetEnumEntries().end());
    return all_scales;
}

const std::set<UncertaintyScale>& GetActiveUncertaintyScales(UncertaintySource unc_source)
{
    static const std::set<UncertaintyScale> none_scales = { UncertaintyScale::Central };
    static const std::set<UncertaintyScale> up_down_scales = { UncertaintyScale::Up, UncertaintyScale::Down };
    return unc_source == UncertaintySource::None ? none_scales : up_down_scales;
}

} // namespace analysis
