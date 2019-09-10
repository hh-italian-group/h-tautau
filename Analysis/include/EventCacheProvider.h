/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Analysis/include/EventCandidate.h"
#include "SVfitAnaInterface.h"
#include "KinFitInterface.h"

namespace analysis {

class EventCacheProvider {
public:
    using LegPair = ntuple::LegPair;

    template<typename Event>
    EventCacheProvider(const Event event);

    bool TryGetKinFit(kin_fit::FitResults kinfit_result, LegPair htt_pair, LegPair hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);
    bool TryGetSVFit(sv_fit_ana::FitResults kinfit_result, LegPair htt_pair, LegPair hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);

    struct EventCacheKey{
        LegPair htt_pair;
        LegPair hbb_pair;
        UncertaintySource unc_source;
        UncertaintyScale unc_scale;

    };

private:
    std::map<EventCacheKey,kin_fit::FitResults> kinFit_map;
    std::map<EventCacheKey,sv_fit_ana::FitResults> SVFit_map;

};

} // namespace analysis
