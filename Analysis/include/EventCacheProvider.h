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

    EventCacheProvider(){} //default constructor

    template<typename Event>
    EventCacheProvider(const Event& event){
        AddEvent(event);
    }

    template<typename Event>
    void AddEvent(const Event& event)
    {
        for(unsigned n = 0; n < event.kinFit_Higgs_index.size(); ++n){
            KinFitKey kinFitKey(ntuple::LegIndexToPair(event.kinFit_Higgs_index.at(n)),
                                ntuple::LegIndexToPair(event.kinFit_jetPairId.at(n)),
                                static_cast<UncertaintySource>(event.kinFit_unc_source.at(n)),
                                static_cast<UncertaintyScale>(event.kinFit_unc_scale.at(n)));
            kin_fit::FitResults kinFit_results(event.kinFit_m.at(n),event.kinFit_chi2.at(n),
                                                  0,event.kinFit_convergence.at(n));
            kinFit_map[kinFitKey] = kinFit_results;
        }

        for(unsigned n = 0; n < event.SVfit_Higgs_index.size(); ++n){
            SVFitKey SVFitKey(ntuple::LegIndexToPair(event.SVfit_Higgs_index.at(n)),
                              static_cast<UncertaintySource>(event.SVfit_unc_source.at(n)),
                              static_cast<UncertaintyScale>(event.SVfit_unc_scale.at(n)));
            sv_fit_ana::FitResults SVFit_results(event.SVfit_is_valid.at(n),
                                                    analysis::LorentzVectorM(event.SVfit_p4.at(n)),
                                                    analysis::LorentzVectorM(event.SVfit_p4_error.at(n)),
                                                    event.SVfit_mt.at(n),
                                                    event.SVfit_mt_error.at(n));
            SVFit_map[SVFitKey] = SVFit_results;
        }
        std::cout << "event.jet_hh_score_index.size(): " << event.jet_hh_score_index.size() << std::endl;
        for(unsigned n = 0; n < event.jet_hh_score_index.size(); ++n){
            JetScoreKey jetScoreKey(event.jet_hh_score_index.at(n),
                              static_cast<UncertaintySource>(event.jet_hh_score_unc_source.at(n)),
                              static_cast<UncertaintyScale>(event.jet_hh_score_unc_scale.at(n)));
            jetScore_map[jetScoreKey] = event.jet_hh_score_value.at(n);
        }
    }

    template<typename Event>
    void FillEvent(Event& event) const
    {
        event.SVfit_Higgs_index.clear();
        event.SVfit_is_valid.clear();
        event.SVfit_p4.clear();
        event.SVfit_p4_error.clear();
        event.SVfit_mt.clear();
        event.SVfit_mt_error.clear();
        event.SVfit_unc_source.clear();
        event.SVfit_unc_scale.clear();
        for(const auto& iter_sv : SVFit_map){
            const sv_fit_ana::FitResults& SVFit_results = iter_sv.second;
            event.SVfit_Higgs_index.push_back(ntuple::LegPairToIndex(iter_sv.first.htt_pair));
            event.SVfit_is_valid.push_back(SVFit_results.has_valid_momentum);
            event.SVfit_p4.push_back(static_cast<ntuple::LorentzVectorM>(SVFit_results.momentum));
            event.SVfit_p4_error.push_back(static_cast<ntuple::LorentzVectorM>(SVFit_results.momentum_error));
            event.SVfit_mt.push_back(static_cast<Float_t>(SVFit_results.transverseMass));
            event.SVfit_mt_error.push_back(static_cast<Float_t>(SVFit_results.transverseMass_error));
            event.SVfit_unc_source.push_back(static_cast<Int_t>(iter_sv.first.unc_source));
            event.SVfit_unc_scale.push_back(static_cast<Int_t>(iter_sv.first.unc_scale));
        }

        event.kinFit_Higgs_index.clear();
        event.kinFit_jetPairId.clear();
        event.kinFit_m.clear();
        event.kinFit_chi2.clear();
        event.kinFit_convergence.clear();
        event.kinFit_unc_source.clear();
        event.kinFit_unc_scale.clear();
        for(const auto& iter_kf : kinFit_map){
            const kin_fit::FitResults& kinFit_results = iter_kf.second;
            event.kinFit_Higgs_index.push_back(static_cast<size_t>(ntuple::LegPairToIndex(iter_kf.first.htt_pair)));
            event.kinFit_jetPairId.push_back(static_cast<size_t>(ntuple::LegPairToIndex(iter_kf.first.hbb_pair)));
            event.kinFit_m.push_back(static_cast<Float_t>(kinFit_results.mass));
            event.kinFit_chi2.push_back(static_cast<Float_t>(kinFit_results.chi2));
            event.kinFit_convergence.push_back(kinFit_results.convergence);
            event.kinFit_unc_source.push_back(static_cast<Int_t>(iter_kf.first.unc_source));
            event.kinFit_unc_scale.push_back(static_cast<Int_t>(iter_kf.first.unc_scale));
        }

        event.jet_hh_score_index.clear();
        event.jet_hh_score_unc_scale.clear();
        event.jet_hh_score_unc_source.clear();
        event.jet_hh_score_value.clear();
        for(const auto& iter_jet : jetScore_map){
            event.jet_hh_score_index.push_back(iter_jet.first.jet_index);
            event.jet_hh_score_unc_scale.push_back(static_cast<Int_t>(iter_jet.first.unc_scale));
            event.jet_hh_score_unc_source.push_back(static_cast<Int_t>(iter_jet.first.unc_source));
            event.jet_hh_score_value.push_back(iter_jet.second);
        }
    }

    bool TryGetKinFit(kin_fit::FitResults& kinfit_result, const LegPair& htt_pair, const LegPair& hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);
    bool TryGetSVFit(sv_fit_ana::FitResults& svfit_result, const LegPair& htt_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);

    struct SVFitKey{
        LegPair htt_pair;
        UncertaintySource unc_source;
        UncertaintyScale unc_scale;

        SVFitKey();
        SVFitKey(LegPair _htt_pair,UncertaintySource _unc_source,UncertaintyScale _unc_scale);

        bool operator<(const SVFitKey& other) const;

        virtual ~SVFitKey(){} // virtual destructor
        SVFitKey( const SVFitKey & ) = default;
    };

    struct KinFitKey : SVFitKey {
        LegPair hbb_pair;

        KinFitKey();
        KinFitKey(LegPair _htt_pair, LegPair _hbb_pair,UncertaintySource _unc_source,UncertaintyScale _unc_scale);
        bool operator<(const KinFitKey& other) const;
    };

    struct JetScoreKey{
        size_t jet_index;
        UncertaintySource unc_source;
        UncertaintyScale unc_scale;

        JetScoreKey();
        JetScoreKey(size_t _jet_index,UncertaintySource _unc_source,UncertaintyScale _unc_scale);

        bool operator<(const JetScoreKey& other) const;

        virtual ~JetScoreKey(){} // virtual destructor
        JetScoreKey( const JetScoreKey & ) = default;
    };



private:
    std::map<KinFitKey,kin_fit::FitResults> kinFit_map;
    std::map<SVFitKey,sv_fit_ana::FitResults> SVFit_map;
    std::map<JetScoreKey,float> jetScore_map;

};

} // namespace analysis
