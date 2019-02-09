/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/TauIdResults.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

namespace tau_id {

TauIdDescriptor::TauIdDescriptor(TauIdDiscriminator _discriminator, const std::string& _name_pattern, const std::string& wp_list)
    : discriminator(_discriminator), name_pattern(_name_pattern)
{
    auto wp_names = SplitValueList(wp_list, false, ", \t", true);
    for(const auto& wp_name : wp_names)
        working_points.push_back(analysis::Parse<DiscriminatorWP>(wp_name));
}

std::string TauIdDescriptor::ToString(DiscriminatorWP wp) const
{
    std::string name = name_pattern;
    boost::algorithm::replace_all(name, "{wp}", analysis::ToString(wp));
    boost::algorithm::replace_all(name, "{raw}", "");
    boost::algorithm::replace_all(name, "{Raw}", "");
    return name;
}

std::string TauIdDescriptor::ToStringRaw() const
{
    std::string name = name_pattern;
    boost::algorithm::replace_all(name, "{wp}", "");
    boost::algorithm::replace_all(name, "{raw}", "raw");
    boost::algorithm::replace_all(name, "{Raw}", "Raw");
    return name;
}

}

TauIdResults::ResultDescriptor::ResultDescriptor(TauIdDiscriminator _discriminator, DiscriminatorWP _wp) :
    discriminator(_discriminator), wp(_wp) {}

bool TauIdResults::ResultDescriptor::operator<(const ResultDescriptor& other) const
{
    if(discriminator != other.discriminator) return discriminator < other.discriminator;
    return wp < other.wp;
}

std::string TauIdResults::ResultDescriptor::ToString() const
{
    return tau_id::GetTauIdDescriptors().at(discriminator).ToString(wp);
}

const TauIdResults::ResultDescriptorCollection& TauIdResults::GetResultDescriptors()
{
    static const auto& discriminators = tau_id::GetOrderedTauIdDiscriminators();
    static const auto& id_descriptors = tau_id::GetTauIdDescriptors();

    auto createDescriptors = [&]() {
        ResultDescriptorCollection descs;
        for(const auto& discriminator : discriminators) {
            for(auto wp : id_descriptors.at(discriminator).working_points)
                descs.emplace_back(discriminator, wp);
        }
        return descs;
    };

    static const ResultDescriptorCollection descriptors = createDescriptors();
    return descriptors;
}

const TauIdResults::BitRefByDescCollection& TauIdResults::GetBitRefsByDesc()
{
    auto createBitRefs = []() {
        BitRefByDescCollection bit_refs;
        const auto& descs = GetResultDescriptors();
        for(size_t n = 0; n < descs.size(); ++n) {
            const auto& desc = descs.at(n);
            if(bit_refs.count(desc))
                throw exception("Duplicated descriptor of tau ID result = '%1%'") % desc.ToString();
            bit_refs[desc] = n;
        }
        return bit_refs;
    };
    static const BitRefByDescCollection bit_refs_by_desc = createBitRefs();
    return bit_refs_by_desc;
}

const TauIdResults::BitRefByNameCollection& TauIdResults::GetBitRefsByName()
{
    auto createBitRefs = []() {
        BitRefByNameCollection bit_refs;
        const auto& descs = GetBitRefsByDesc();
        for(const auto& item : descs)
            bit_refs[item.first.ToString()] = item.second;
        return bit_refs;
    };
    static const BitRefByNameCollection bit_refs_by_name = createBitRefs();
    return bit_refs_by_name;
}

TauIdResults::TauIdResults() : result_bits(0) {}
TauIdResults::TauIdResults(BitsContainer _result_bits) : result_bits(_result_bits) {}

TauIdResults::BitsContainer TauIdResults::GetResultBits() const { return result_bits.to_ullong(); }

bool TauIdResults::Result(size_t index) const { CheckIndex(index); return result_bits[index]; }
void TauIdResults::SetResult(size_t index, bool value) { CheckIndex(index); result_bits[index] = value; }

bool TauIdResults::Result(TauIdDiscriminator discriminator, DiscriminatorWP wp) const
{
    const ResultDescriptor desc(discriminator, wp);
    const auto& bit_refs = GetBitRefsByDesc();
    auto iter = bit_refs.find(desc);
    if(iter == bit_refs.end())
        throw exception("Result bit not found for %1% working point of %2%.") % wp % discriminator;
    return Result(iter->second);
}

bool TauIdResults::Result(const std::string& name) const
{
    const auto& bit_refs = GetBitRefsByName();
    auto iter = bit_refs.find(name);
    if(iter == bit_refs.end())
        throw exception("Result bit not found for '%1%'.") % name;
    return Result(iter->second);
}

void TauIdResults::CheckIndex(size_t index) const
{
    if(index >= MaxNumberOfIds)
        throw exception("Tau ID index is out of range.");
}

} // namespace analysis
