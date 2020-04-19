/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "h-tautau/Analysis/include/SVfitAnaInterface.h"
#include "h-tautau/Core/include/CacheTuple.h"

namespace analysis {

class EventCacheProvider {
public:
    using LegPair = ntuple::LegPair;

    struct SVFitKey {
        size_t htt_index;
        UncertaintySource unc_source;
        UncertaintyScale unc_scale;

        SVFitKey(){}
        SVFitKey(const SVFitKey&) = default;
        SVFitKey(size_t _htt_index, UncertaintySource _unc_source, UncertaintyScale _unc_scale);
        bool operator<(const SVFitKey& other) const;
        virtual ~SVFitKey(){}
    };

    struct KinFitKey : SVFitKey {
        size_t hbb_index;

        KinFitKey(){}
        KinFitKey(size_t _htt_index, size_t _hbb_index, UncertaintySource _unc_source, UncertaintyScale _unc_scale);
        bool operator<(const KinFitKey& other) const;
    };

    struct HHBtagKey : SVFitKey {
        size_t jet_index;

        HHBtagKey(){}
        HHBtagKey(size_t _htt_index, size_t _jet_index, UncertaintySource _unc_source, UncertaintyScale _unc_scale);
        bool operator<(const HHBtagKey& other) const;
    };

    EventCacheProvider() {}
    template<typename Event>
    EventCacheProvider(const Event& event) { AddEvent(event); }

    void AddSVfitResults(size_t htt_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                         const sv_fit_ana::FitResults& fit_results);
    void AddKinFitResults(size_t htt_index, size_t hbb_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                          const kin_fit::FitResults& fit_results);
    void AddHHbtagResults(size_t htt_index, size_t jet_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                          float hhbtag_score);
    void AddHHbtagResults(size_t htt_index, UncertaintySource unc_source, UncertaintyScale unc_scale,
                          const std::vector<float>& hhbtag_scores);

    template<typename Event>
    void AddEvent(const Event& event, bool add_svfit = true, bool add_kinfit = true, bool add_hhbtag = true)
    {
        if(add_svfit) {
            for(size_t n = 0; n < event.SVfit_htt_index.size(); ++n) {
                const sv_fit_ana::FitResults fit_results(event.SVfit_is_valid.at(n),
                                                         LorentzVectorM(event.SVfit_p4.at(n)),
                                                         LorentzVectorM(event.SVfit_p4_error.at(n)),
                                                         event.SVfit_mt.at(n), event.SVfit_mt_error.at(n));
                AddSVfitResults(event.SVfit_htt_index.at(n),
                                static_cast<UncertaintySource>(event.SVfit_unc_source.at(n)),
                                static_cast<UncertaintyScale>(event.SVfit_unc_scale.at(n)), fit_results);
            }
        }

        if(add_kinfit) {
            for(size_t n = 0; n < event.kinFit_htt_index.size(); ++n) {
                const kin_fit::FitResults fit_results(event.kinFit_m.at(n),event.kinFit_chi2.at(n), 0,
                                                      event.kinFit_convergence.at(n));
                AddKinFitResults(event.kinFit_htt_index.at(n), event.kinFit_hbb_index.at(n),
                                 static_cast<UncertaintySource>(event.kinFit_unc_source.at(n)),
                                 static_cast<UncertaintyScale>(event.kinFit_unc_scale.at(n)), fit_results);
            }
        }

        if(add_hhbtag) {
            for(unsigned n = 0; n < event.jet_HHbtag_htt_index.size(); ++n) {
                AddHHbtagResults(event.jet_HHbtag_htt_index.at(n), event.jet_HHbtag_jet_index.at(n),
                                 static_cast<UncertaintySource>(event.jet_HHbtag_unc_source.at(n)),
                                 static_cast<UncertaintyScale>(event.jet_HHbtag_unc_scale.at(n)),
                                 event.jet_HHbtag_value.at(n));
            }
        }
    }

    template<typename Event>
    void FillEvent(Event& event, bool fill_svfit = true, bool fill_kinfit = true, bool fill_hhbtag = true) const
    {
        using LorentzVectorM = typename decltype(event.SVfit_p4)::value_type;

        if(fill_svfit) {
            event.SVfit_htt_index.clear();
            event.SVfit_is_valid.clear();
            event.SVfit_p4.clear();
            event.SVfit_p4_error.clear();
            event.SVfit_mt.clear();
            event.SVfit_mt_error.clear();
            event.SVfit_unc_source.clear();
            event.SVfit_unc_scale.clear();
            for(const auto& [key, result] : SVFit_map){
                event.SVfit_htt_index.push_back(static_cast<UInt_t>(key.htt_index));
                event.SVfit_is_valid.push_back(result.has_valid_momentum);
                event.SVfit_p4.push_back(LorentzVectorM(result.momentum));
                event.SVfit_p4_error.push_back(LorentzVectorM(result.momentum_error));
                event.SVfit_mt.push_back(static_cast<Float_t>(result.transverseMass));
                event.SVfit_mt_error.push_back(static_cast<Float_t>(result.transverseMass_error));
                event.SVfit_unc_source.push_back(static_cast<Int_t>(key.unc_source));
                event.SVfit_unc_scale.push_back(static_cast<Int_t>(key.unc_scale));
            }
        }

        if(fill_kinfit) {
            event.kinFit_htt_index.clear();
            event.kinFit_hbb_index.clear();
            event.kinFit_unc_source.clear();
            event.kinFit_unc_scale.clear();
            event.kinFit_m.clear();
            event.kinFit_chi2.clear();
            event.kinFit_convergence.clear();
            for(const auto& [key, result] : kinFit_map){
                event.kinFit_htt_index.push_back(static_cast<UInt_t>(key.htt_index));
                event.kinFit_hbb_index.push_back(static_cast<UInt_t>(key.hbb_index));
                event.kinFit_unc_source.push_back(static_cast<Int_t>(key.unc_source));
                event.kinFit_unc_scale.push_back(static_cast<Int_t>(key.unc_scale));
                event.kinFit_m.push_back(static_cast<Float_t>(result.mass));
                event.kinFit_chi2.push_back(static_cast<Float_t>(result.chi2));
                event.kinFit_convergence.push_back(result.convergence);
            }
        }

        if(fill_hhbtag) {
            event.jet_HHbtag_htt_index.clear();
            event.jet_HHbtag_jet_index.clear();
            event.jet_HHbtag_unc_source.clear();
            event.jet_HHbtag_unc_scale.clear();
            event.jet_HHbtag_value.clear();
            for(const auto& [key, result] : hhBtag_map) {
                event.jet_HHbtag_htt_index.push_back(static_cast<UInt_t>(key.htt_index));
                event.jet_HHbtag_jet_index.push_back(static_cast<UInt_t>(key.jet_index));
                event.jet_HHbtag_unc_source.push_back(static_cast<Int_t>(key.unc_scale));
                event.jet_HHbtag_unc_scale.push_back(static_cast<Int_t>(key.unc_source));
                event.jet_HHbtag_value.push_back(result);
            }
        }
    }

    bool IsEmpty() const;

    boost::optional<sv_fit_ana::FitResults> TryGetSVFit(size_t htt_index, UncertaintySource unc_source,
                                                        UncertaintyScale unc_scale) const;
    boost::optional<kin_fit::FitResults> TryGetKinFit(size_t htt_index, size_t hbb_index, UncertaintySource unc_source,
                                                      UncertaintyScale unc_scale) const;
    boost::optional<float> TryGetHHbtag(size_t htt_index, size_t jet_index, UncertaintySource unc_source,
                                        UncertaintyScale unc_scale) const;

private:
    std::map<SVFitKey, sv_fit_ana::FitResults> SVFit_map;
    std::map<KinFitKey, kin_fit::FitResults> kinFit_map;
    std::map<HHBtagKey, float> hhBtag_map;
};

class EventCacheSource {
public:
    EventCacheSource(const std::string& file_name, const std::string& tree_name);
    bool Read(Long64_t entry_index, EventCacheProvider& provider);
    boost::optional<Long64_t> GetCurrentEntryIndex() const;
    size_t GetTotalNumberOfEntries() const;
    size_t GetRemainingNumberOfEntries() const;
    const std::vector<cache_tuple::CacheProdSummary>& GetSummary() const;

private:
    void Advance();

private:
    std::shared_ptr<TFile> file;
    std::shared_ptr<cache_tuple::CacheTuple> cache;
    Long64_t current_cache_entry;
    std::vector<cache_tuple::CacheProdSummary> summary;
};

class EventCacheReader {
public:
    EventCacheReader(const std::vector<std::string>& cache_files, const std::string& tree_name);
    EventCacheProvider Read(Long64_t entry_index);
    boost::optional<Long64_t> GetCurrentEntryIndex() const;
    size_t GetTotalNumberOfEntries() const;
    size_t GetRemainingNumberOfEntries() const;
    const std::vector<cache_tuple::CacheProdSummary>& GetSummary() const;

private:
    std::vector<EventCacheSource> sources;
    Long64_t last_entry_index;
    size_t n_total, n_remaining;
    std::vector<cache_tuple::CacheProdSummary> summary;
};

} // namespace analysis
