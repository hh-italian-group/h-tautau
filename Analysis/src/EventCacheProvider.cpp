/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventCacheProvider.h"

namespace analysis {

    bool EventCacheProvider::TryGetKinFit(kin_fit::FitResults kinfit_result, LegPair htt_pair,
            LegPair hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale)
    {
        KinFitKey kinFitKey(htt_pair,hbb_pair,unc_source,unc_scale);
        bool gotKinFit = false;
        if(kinFit_map.count(kinFitKey)){
            gotKinFit = true;
            kinfit_result = kinFit_map.at(kinFitKey);
        }
        return gotKinFit;
    }

    bool EventCacheProvider::TryGetSVFit(sv_fit_ana::FitResults svfit_result, LegPair htt_pair,
        UncertaintySource unc_source, UncertaintyScale unc_scale)
    {
        SVFitKey svFitKey(htt_pair,unc_source,unc_scale);
        bool gotSVFit = false;
        if(SVFit_map.count(svFitKey)){
            gotSVFit = true;
            svfit_result = SVFit_map.at(svFitKey);
        }
        return gotSVFit;
    }


} // namespace analysis
