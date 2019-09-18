/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventCacheProvider.h"

namespace analysis {

    bool EventCacheProvider::TryGetKinFit(kin_fit::FitResults& kinfit_result, const LegPair& htt_pair,
        const LegPair& hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale)
    {
        KinFitKey kinFitKey(htt_pair,hbb_pair,unc_source,unc_scale);
        bool gotKinFit = false;
        auto iter = kinFit_map.find(kinFitKey);
        if(iter != kinFit_map.end()){
            gotKinFit = true;
            kinfit_result = iter->second;
        }
        return gotKinFit;
    }

    bool EventCacheProvider::TryGetSVFit(sv_fit_ana::FitResults& svfit_result, const LegPair& htt_pair,
        UncertaintySource unc_source, UncertaintyScale unc_scale)
    {
        SVFitKey svFitKey(htt_pair,unc_source,unc_scale);
        bool gotSVFit = false;
        auto iter = SVFit_map.find(svFitKey);
        if(iter != SVFit_map.end()){
            gotSVFit = true;
            svfit_result = iter->second;
        }
        return gotSVFit;
    }

    EventCacheProvider::SVFitKey::SVFitKey() : htt_pair(ntuple::UndefinedLegPair()),
            unc_source(UncertaintySource::None), unc_scale(UncertaintyScale::Central) { }

    EventCacheProvider::SVFitKey::SVFitKey(LegPair _htt_pair,UncertaintySource _unc_source,UncertaintyScale _unc_scale) :
            htt_pair(_htt_pair), unc_source(_unc_source), unc_scale(_unc_scale) { }

    bool EventCacheProvider::SVFitKey::operator<(const SVFitKey& other) const
    {
        if(htt_pair != other.htt_pair) return htt_pair < other.htt_pair;
        if(unc_source != other.unc_source) return unc_source < other.unc_source;
        return unc_scale < other.unc_scale;
    }

    EventCacheProvider::KinFitKey::KinFitKey() : SVFitKey(),hbb_pair(ntuple::UndefinedLegPair()) { }

    EventCacheProvider::KinFitKey::KinFitKey(LegPair _htt_pair, LegPair _hbb_pair, UncertaintySource _unc_source,UncertaintyScale _unc_scale) :
            SVFitKey(_htt_pair,_unc_source,_unc_scale), hbb_pair(_hbb_pair) { }

    bool EventCacheProvider::KinFitKey::operator<(const KinFitKey& other) const
    {
        if(htt_pair != other.htt_pair) return htt_pair < other.htt_pair;
        if(hbb_pair != other.hbb_pair) return hbb_pair < other.hbb_pair;
        if(unc_source != other.unc_source) return unc_source < other.unc_source;
        return unc_scale < other.unc_scale;
    }

} // namespace analysis
