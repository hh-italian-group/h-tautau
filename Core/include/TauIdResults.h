/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include "AnalysisTypes.h"

namespace analysis {

#define TAU_IDS() \
    TAU_ID(againstElectronMVA6, "againstElectron{wp}MVA6{Raw}", "VLoose Loose Medium Tight VTight") \
    TAU_ID(againstMuon3, "againstMuon{wp}3", "Loose Tight") \
    TAU_ID(byCombinedIsolationDeltaBetaCorr3Hits, "by{wp}CombinedIsolationDeltaBetaCorr{Raw}3Hits", \
           "Loose Medium Tight") \
    TAU_ID(byPhotonPtSumOutsideSignalCone, "byPhotonPtSumOutsideSignalCone", "Medium") \
    TAU_ID(byIsolationMVArun2v1DBoldDMwLT, "by{wp}IsolationMVArun2v1DBoldDMwLT{raw}", \
           "VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2v1DBdR03oldDMwLT, "by{wp}IsolationMVArun2v1DBdR03oldDMwLT{raw}", \
           "VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2v1DBoldDMwLT2016, "by{wp}IsolationMVArun2v1DBoldDMwLT{raw}2016", \
           "VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBoldDMwLT2017, "by{wp}IsolationMVArun2017v2DBoldDMwLT{raw}2017", \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    TAU_ID(byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, "by{wp}IsolationMVArun2017v2DBoldDMdR0p3wLT{raw}2017", \
           "VVLoose VLoose Loose Medium Tight VTight VVTight") \
    /**/

#define TAU_ID(name, pattern, wp_list) name,
enum class TauIdDiscriminator { TAU_IDS() };
#undef TAU_ID

#define TAU_ID(name, pattern, wp_list) TauIdDiscriminator::name,
namespace tau_id {
inline const std::vector<TauIdDiscriminator>& GetOrderedTauIdDiscriminators()
{
    static const std::vector<TauIdDiscriminator> ordered_tau_ids = { TAU_IDS() };
    return ordered_tau_ids;
}
}
#undef TAU_ID

#define TAU_ID(name, pattern, wp_list) { TauIdDiscriminator::name, #name },
ENUM_NAMES(TauIdDiscriminator) = { TAU_IDS() };
#undef TAU_ID

namespace tau_id {
struct TauIdDescriptor {
    TauIdDiscriminator discriminator;
    std::string name_pattern;
    std::vector<DiscriminatorWP> working_points;

    TauIdDescriptor(TauIdDiscriminator _discriminator, const std::string& _name_pattern, const std::string& wp_list);
    std::string ToString(DiscriminatorWP wp) const;
    std::string ToStringRaw() const;
};

using TauIdDescriptorCollection = std::map<TauIdDiscriminator, TauIdDescriptor>;

#define TAU_ID(name, pattern, wp_list) \
    { TauIdDiscriminator::name, TauIdDescriptor(TauIdDiscriminator::name, pattern, wp_list) },
inline const TauIdDescriptorCollection& GetTauIdDescriptors()
{
    static const TauIdDescriptorCollection descriptors = { TAU_IDS() };
    return descriptors;
}
#undef TAU_ID

}

#undef TAU_IDS

class TauIdResults {
public:
    using BitsContainer = unsigned long long;
    static constexpr size_t MaxNumberOfIds = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfIds>;

    struct ResultDescriptor {
        TauIdDiscriminator discriminator;
        DiscriminatorWP wp;

        ResultDescriptor() {}
        ResultDescriptor(TauIdDiscriminator _discriminator, DiscriminatorWP _wp);
        bool operator<(const ResultDescriptor& other) const;
        std::string ToString() const;
    };

    using ResultDescriptorCollection = std::vector<ResultDescriptor>;
    using BitRefByDescCollection = std::map<ResultDescriptor, size_t>;
    using BitRefByNameCollection = std::map<std::string, size_t>;

    static const ResultDescriptorCollection& GetResultDescriptors();
    static const BitRefByDescCollection& GetBitRefsByDesc();
    static const BitRefByNameCollection& GetBitRefsByName();

    TauIdResults();
    TauIdResults(BitsContainer _result_bits);

    BitsContainer GetResultBits() const;

    bool Result(size_t index) const;
    void SetResult(size_t index, bool value);
    bool Result(TauIdDiscriminator discriminator, DiscriminatorWP wp) const;
    bool Result(const std::string& name) const;

private:
    void CheckIndex(size_t index) const;

private:
    Bits result_bits;
};

} // namespace analysis
