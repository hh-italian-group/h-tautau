/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include "AnalysisTypes.h"
#include <iostream>
#include "AnalysisTools/Core/include/TextIO.h"
#include "DiscriminatorIdResults.h"

namespace analysis {

#define TAU_IDS() \
    TAU_ID(againstElectronMVA6, "againstElectron{wp}MVA6{Raw}", true, "VLoose Loose Medium Tight VTight") \
    TAU_ID(againstElectronMVA62018, "againstElectron{wp}MVA6{Raw}2018", true, "VLoose Loose Medium Tight VTight") \
    TAU_ID(againstMuon3, "againstMuon{wp}3", false, "Loose Tight") \
    TAU_ID(byCombinedIsolationDeltaBetaCorr3Hits, "by{wp}CombinedIsolationDeltaBetaCorr{Raw}3Hits", true, \
           "Loose Medium Tight") \
    TAU_ID(byIsolationMVArun2v1DBoldDMwLT2016, "by{wp}IsolationMVArun2v1DBoldDMwLT{raw}2016", true, \
           "VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBoldDMwLT2017, "by{wp}IsolationMVArun2017v2DBoldDMwLT{raw}2017", true, \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDeepTau2017v2p1VSe, "by{wp}DeepTau2017v2p1VSe{raw}", true, \
           "VVVLoose VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDeepTau2017v2p1VSmu, "by{wp}DeepTau2017v2p1VSmu{raw}", true, \
           "VLoose Loose Medium Tight") \
    TAU_ID(byDeepTau2017v2p1VSjet, "by{wp}DeepTau2017v2p1VSjet{raw}", true, \
           "VVVLoose VVLoose VLoose Loose Medium Tight VTight VVTight") \
    /**/


#define TAU_ID(name, pattern, has_raw, wp_list) name,
enum class TauIdDiscriminator { TAU_IDS() };
#undef TAU_ID

#define TAU_ID(name, pattern, has_raw, wp_list) TauIdDiscriminator::name,
namespace tau_id {
inline const std::vector<TauIdDiscriminator>& GetOrderedTauIdDiscriminators()
{
    static const std::vector<TauIdDiscriminator> ordered_tau_ids = { TAU_IDS() };
    return ordered_tau_ids;
}
}
#undef TAU_ID

#define TAU_ID(name, pattern, has_raw, wp_list) { TauIdDiscriminator::name, #name },
ENUM_NAMES(TauIdDiscriminator) = { TAU_IDS() };
#undef TAU_ID

namespace tau_id {
struct TauIdDescriptor {
    TauIdDiscriminator discriminator;
    std::string name_pattern;
    bool has_raw;
    std::string raw_name;
    std::map<DiscriminatorWP, std::string> working_points;

    TauIdDescriptor(TauIdDiscriminator _discriminator, const std::string& _name_pattern, bool _has_raw,
                    const std::string& wp_list);

    std::string ToString(DiscriminatorWP wp) const;
    std::string ToStringRaw() const;

    template<typename Tuple, typename Tau>
    void FillTuple(Tuple& tuple, const Tau* tau, float default_value, const std::string& prefix = "",
                   const std::string& raw_suffix = "raw") const
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        if(has_raw){
            float value =  tau && tau->isTauIDAvailable(raw_name) ? tau->tauID(raw_name) : default_value;
            tuple.template get<std::vector<float_t>>(prefix + disc_name + raw_suffix).push_back(value);
        }
        if(!working_points.empty()) {
            DiscriminatorIdResults id_results;
            for(const auto& wp_entry : working_points) {
                const bool result = tau && tau->isTauIDAvailable(wp_entry.second) && tau->tauID(wp_entry.second) > 0.5;
                id_results.SetResult(wp_entry.first, result);
            }
            tuple.template get<std::vector<DiscriminatorIdResults::BitsContainer>>(prefix + disc_name).push_back(id_results.GetResultBits());
        }
    }

    template<typename Tuple>
    DiscriminatorIdResults GetIdResult(Tuple& tuple, size_t n, const std::string& prefix = "") const
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        return DiscriminatorIdResults(tuple.template get<DiscriminatorIdResults::BitsContainer>(prefix + disc_name)->at(n));
    }

    template<typename Tuple>
    float_t GetRawId(Tuple& tuple, size_t n, const std::string& prefix = "", const std::string& raw_suffix = "raw") const
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        return tuple.template get<float>(prefix + disc_name + raw_suffix)->at(n);
    }
};

using TauIdDescriptorCollection = std::map<TauIdDiscriminator, TauIdDescriptor>;

#define TAU_ID(name, pattern, has_raw, wp_list) \
    { TauIdDiscriminator::name, TauIdDescriptor(TauIdDiscriminator::name, pattern, has_raw, wp_list) },
inline const TauIdDescriptorCollection& GetTauIdDescriptors()
{
    static const TauIdDescriptorCollection descriptors = { TAU_IDS() };
    return descriptors;
}
#undef TAU_ID

} // namespace tau_id

} // namespace analysis
