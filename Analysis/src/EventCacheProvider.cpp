/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventCacheProvider.h"

namespace analysis {

EventCacheProvider::SVFitKey::SVFitKey(size_t _htt_index, UncertaintySource _unc_source, UncertaintyScale _unc_scale) :
        htt_index(_htt_index), unc_source(_unc_source), unc_scale(_unc_scale)
{
}

bool EventCacheProvider::SVFitKey::operator<(const SVFitKey& other) const
{
    if(htt_index != other.htt_index) return htt_index < other.htt_index;
    if(unc_source != other.unc_source) return unc_source < other.unc_source;
    return unc_scale < other.unc_scale;
}

EventCacheProvider::KinFitKey::KinFitKey(size_t _htt_index, size_t _hbb_index, UncertaintySource _unc_source,
                                         UncertaintyScale _unc_scale) :
        SVFitKey(_htt_index, _unc_source, _unc_scale), hbb_index(_hbb_index)
{
}

bool EventCacheProvider::KinFitKey::operator<(const KinFitKey& other) const
{
    if(hbb_index != other.hbb_index) return hbb_index < other.hbb_index;
    return SVFitKey::operator<(other);
}

EventCacheProvider::HHBtagKey::HHBtagKey(size_t _htt_index, size_t _jet_index, UncertaintySource _unc_source,
                                         UncertaintyScale _unc_scale) :
        SVFitKey(_htt_index, _unc_source, _unc_scale), jet_index(_jet_index)
{
}

bool EventCacheProvider::HHBtagKey::HHBtagKey::operator<(const HHBtagKey& other) const
{
    if(jet_index != other.jet_index) return jet_index < other.jet_index;
    return SVFitKey::operator<(other);
}


void EventCacheProvider::AddSVfitResults(size_t htt_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                                         const sv_fit_ana::FitResults& fit_results)
{
    const SVFitKey key(htt_index, unc_source, unc_scale);
    if(SVFit_map.count(key))
        throw exception("EventCacheProvider: duplicated SVfit entry.");
    SVFit_map[key] = fit_results;
}

void EventCacheProvider::AddKinFitResults(size_t htt_index, size_t hbb_index, UncertaintySource unc_source,
                                          UncertaintyScale unc_scale, const kin_fit::FitResults& fit_results)
{
    const KinFitKey key(htt_index, hbb_index, unc_source, unc_scale);
    if(kinFit_map.count(key))
        throw exception("EventCacheProvider: duplicated KinFit entry.");
    kinFit_map[key] = fit_results;
}

void EventCacheProvider::AddHHbtagResults(size_t htt_index, size_t jet_index, UncertaintySource unc_source,
                                          UncertaintyScale unc_scale, float hhbtag_score)
{
    const HHBtagKey key(htt_index, jet_index, unc_source, unc_scale);
    if(hhBtag_map.count(key))
        throw exception("EventCacheProvider: duplicated HHBtag entry.");
    hhBtag_map[key] = hhbtag_score;
}

void EventCacheProvider::AddHHbtagResults(size_t htt_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                                          const std::vector<float>& hhbtag_scores)
{
    for(size_t jet_index = 0; jet_index < hhbtag_scores.size(); ++jet_index)
        AddHHbtagResults(htt_index, jet_index, unc_source, unc_scale, hhbtag_scores.at(jet_index));
}

bool EventCacheProvider::IsEmpty() const
{
    return SVFit_map.empty() && kinFit_map.empty() && hhBtag_map.empty();
}

boost::optional<sv_fit_ana::FitResults> EventCacheProvider::TryGetSVFit(size_t htt_index, UncertaintySource unc_source,
                                                                        UncertaintyScale unc_scale) const
{
    boost::optional<sv_fit_ana::FitResults> result;
    const SVFitKey key(htt_index, unc_source, unc_scale);
    auto iter = SVFit_map.find(key);
    if(iter != SVFit_map.end())
        result = iter->second;
    return result;
}

boost::optional<kin_fit::FitResults> EventCacheProvider::TryGetKinFit(size_t htt_index, size_t hbb_index,
                                                                      UncertaintySource unc_source,
                                                                      UncertaintyScale unc_scale) const
{
    boost::optional<kin_fit::FitResults> result;
    const KinFitKey key(htt_index, hbb_index, unc_source, unc_scale);
    auto iter = kinFit_map.find(key);
    if(iter != kinFit_map.end())
        result = iter->second;
    return result;
}

boost::optional<float> EventCacheProvider::TryGetHHbtag(size_t htt_index, size_t jet_index,
                                                        UncertaintySource unc_source,
                                                        UncertaintyScale unc_scale) const
{
    boost::optional<float> result;
    const HHBtagKey key(htt_index, jet_index, unc_source, unc_scale);
    auto iter = hhBtag_map.find(key);
    if(iter != hhBtag_map.end())
        result = iter->second;
    return result;
}

EventCacheSource::EventCacheSource(const std::string& file_name, const std::string& tree_name) :
    file(root_ext::OpenRootFile(file_name)),
    cache(std::make_shared<cache_tuple::CacheTuple>(tree_name, file.get(), true)),
    current_cache_entry(-1)
{
    Advance();
    cache_tuple::CacheSummaryTuple summary_tuple("summary", file.get(), true);
    for(const auto& s : summary_tuple)
        summary.push_back(s);
}

bool EventCacheSource::Read(Long64_t entry_index, EventCacheProvider& provider)
{
    while(current_cache_entry < cache->GetEntries() && cache->data().entry_index < entry_index)
        Advance();
    if(current_cache_entry < cache->GetEntries() && cache->data().entry_index == entry_index) {
        provider.AddEvent(cache->data());
        Advance();
        return true;
    }
    return false;
}

boost::optional<Long64_t> EventCacheSource::GetCurrentEntryIndex() const
{
    boost::optional<Long64_t> entry_index;
    if(current_cache_entry < cache->GetEntries())
        entry_index = cache->data().entry_index;
    return entry_index;
}

size_t EventCacheSource::GetTotalNumberOfEntries() const { return static_cast<size_t>(cache->GetEntries()); }
size_t EventCacheSource::GetRemainingNumberOfEntries() const
{
    const Long64_t delta = std::max<Long64_t>(cache->GetEntries() - current_cache_entry, 0);
    return static_cast<size_t>(delta);
}
const std::vector<cache_tuple::CacheProdSummary>& EventCacheSource::GetSummary() const { return summary; }

void EventCacheSource::Advance()
{
    ++current_cache_entry;
    if(current_cache_entry < cache->GetEntries())
        cache->GetEntry(current_cache_entry);
}

EventCacheReader::EventCacheReader(const std::vector<std::string>& cache_files, const std::string& tree_name)
    : last_entry_index(-1), n_total(0)
{
    for(const std::string& file_name : cache_files) {
        sources.emplace_back(file_name, tree_name);
        n_total += sources.back().GetTotalNumberOfEntries();
        summary.insert(summary.end(), sources.back().GetSummary().begin(), sources.back().GetSummary().end());
    }
    n_remaining = n_total;
}

EventCacheProvider EventCacheReader::Read(Long64_t entry_index)
{
    if(entry_index <= last_entry_index)
        throw exception("EventCacheReader: only the forward reading is supported.");
    last_entry_index = entry_index;
    EventCacheProvider provider;
    n_remaining = 0;
    for(auto& source : sources) {
        source.Read(entry_index, provider);
        n_remaining += source.GetRemainingNumberOfEntries();
    }
    return provider;
}

boost::optional<Long64_t> EventCacheReader::GetCurrentEntryIndex() const
{
    boost::optional<Long64_t> entry_index;
    for(const auto& source : sources) {
        const auto source_entry_index = source.GetCurrentEntryIndex();
        if(source_entry_index && (!entry_index || *entry_index < *source_entry_index))
            entry_index = source_entry_index;
    }
    return entry_index;
}

size_t EventCacheReader::GetTotalNumberOfEntries() const { return n_total; }
size_t EventCacheReader::GetRemainingNumberOfEntries() const { return n_remaining; }
const std::vector<cache_tuple::CacheProdSummary>& EventCacheReader::GetSummary() const { return summary; }

} // namespace analysis
