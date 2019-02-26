/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include "AnalysisTypes.h"
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
    TAU_ID(byIsolationMVArun2v1DBnewDMwLT2016, "by{wp}IsolationMVArun2v1DBnewDMwLT{raw}2016", true, \
           "VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBoldDMwLT2017, "by{wp}IsolationMVArun2017v2DBoldDMwLT{raw}2017", true, \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, "by{wp}IsolationMVArun2017v2DBoldDMdR0p3wLT{raw}2017", true, \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBnewDMwLT2017, "by{wp}IsolationMVArun2017v2DBnewDMwLT{raw}2017", true, \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDeepTau2017v1VSe, "by{wp}DeepTau2017v1VSe{raw}", true, \
           "VVVLoose VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDeepTau2017v1VSmu, "by{wp}DeepTau2017v1VSmu{raw}", true, \
           "VVVLoose VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDeepTau2017v1VSjet, "by{wp}DeepTau2017v1VSjet{raw}", true, \
           "VVVLoose VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byDpfTau2016v0VSall, "by{wp}DpfTau2016v0VSall{raw}", true, "Tight") \
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
                    const std::string& wp_list) :
        discriminator(_discriminator), name_pattern(_name_pattern), has_raw(_has_raw)
    {
        if(has_raw)
            raw_name = ToStringRaw();
        auto wp_names = SplitValueList(wp_list, false, ", \t", true);
        for(const auto& wp_name : wp_names) {
            const DiscriminatorWP wp = ::analysis::Parse<DiscriminatorWP>(wp_name);
            working_points[wp] = ToString(wp);
        }
    }

    std::string ToString(DiscriminatorWP wp) const
    {
        std::string name = name_pattern;
        boost::algorithm::replace_all(name, "{wp}", analysis::ToString(wp));
        boost::algorithm::replace_all(name, "{raw}", "");
        boost::algorithm::replace_all(name, "{Raw}", "");
        return name;
    }

    std::string ToStringRaw() const
    {
        std::string name = name_pattern;
        boost::algorithm::replace_all(name, "{wp}", "");
        boost::algorithm::replace_all(name, "{raw}", "raw");
        boost::algorithm::replace_all(name, "{Raw}", "Raw");
        return name;
    }

    template<typename Tuple, typename Tau>
    void FillTuple(Tuple& tuple, const Tau* tau, float default_value, const std::string& prefix = "",
                   const std::string& raw_suffix = "raw") const
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        if(has_raw){
            float value =  tau ? tau->tauID(raw_name) : default_value;
            tuple.template get<float>(prefix + disc_name + raw_suffix).push_back(value);
        }
        if(!working_points.empty()) {
            DiscriminatorIdResults id_results;
            for(const auto& wp_entry : working_points) {
                const bool result = tau && tau->tauID(wp_entry.second) > 0.5;
                id_results.SetResult(wp_entry.first, result);
            }
            tuple.template get<DiscriminatorIdResults::BitsContainer>(prefix + disc_name).push_back(id_results.GetResultBits());
        }
    }

    template<typename Tuple>
    DiscriminatorIdResults GetIdResult(Tuple& tuple, int_t n, const std::string& prefix = "")
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        return DiscriminatorIdResults(tuple.template get<DiscriminatorIdResults::BitsContainer>(prefix + disc_name).at(n));
    }

    template<typename Tuple>
    Float_t GetRawId(Tuple& tuple, int_t n, const std::string& prefix = "", const std::string& raw_suffix = "raw")
    {
        const std::string disc_name = ::analysis::ToString(discriminator);
        return tuple.template get<float>(prefix + disc_name + raw_suffix).at(n);
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
