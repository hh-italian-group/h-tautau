/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"
#include "DiscriminatorIdResults.h"
#include "TauIdResults.h"
#include "TriggerResults.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using DiscriminatorIdResults = analysis::DiscriminatorIdResults;
    using DiscriminatorResult = float;
    using Integer = int;
    using RealNumber = float;

    TupleObject(const ntuple::Event& _event);

protected:
    const Event* event;

    static void CheckIndexRange(size_t index, size_t size, std::string_view obj_name, std::string_view branch_name);

    template<typename Value>
    static Value CheckAndGet(size_t index, const std::vector<Value>& col, std::string_view obj_name,
                             std::string_view branch_name)
    {
        CheckIndexRange(index, col.size(), obj_name, branch_name);
        return col.at(index);
    }

    template<typename Value>
    static const Value& CheckAndGetRef(size_t index, const std::vector<Value>& col, std::string_view obj_name,
                                       std::string_view branch_name)
    {
        CheckIndexRange(index, col.size(), obj_name, branch_name);
        return col.at(index);
    }

};

class TupleLepton : public TupleObject {
public:
    TupleLepton(const ntuple::Event& _event, size_t _object_id);
    const LorentzVectorM& p4() const;
    Integer charge() const;
    RealNumber dxy() const;
    RealNumber dz() const;
    RealNumber iso() const;
    analysis::GenLeptonMatch gen_match() const;
    const LorentzVectorM& gen_p4() const;
    Integer decayMode() const;
    analysis::LegType leg_type() const;
    bool passConversionVeto() const;
    bool passEleIsoId(DiscriminatorWP wp) const;
    bool passEleNoIsoId(DiscriminatorWP wp) const;
    bool passMuonId(DiscriminatorWP wp) const;

    bool Passed(analysis::TauIdDiscriminator tauIdDiscriminator, DiscriminatorWP wp) const;
    bool PassedOldDecayMode() const;
    bool PassedNewDecayMode() const;
    DiscriminatorResult GetRawValue(analysis::TauIdDiscriminator tauIdDiscriminator) const;

    int CompareIsolations(const TupleLepton& other, analysis::TauIdDiscriminator disc) const;

private:
    template<typename Value>
    Value CheckAndGet(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGet(object_id, col, "lepton", branch_name);
    }

    template<typename Value>
    const Value& CheckAndGetRef(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGetRef(object_id, col, "lepton", branch_name);
    }

private:
    size_t object_id;
};


class TupleJet : public TupleObject {
public:
    using FilterBits = analysis::TriggerDescriptorCollection::BitsContainer;

    TupleJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    bool PassPuId(DiscriminatorWP wp) const;
    analysis::DiscriminatorIdResults GetPuId() const;
    Float_t GetPuIdRaw() const;
    DiscriminatorResult csv() const;
    DiscriminatorResult deepcsv() const;
    DiscriminatorResult deepFlavour() const;
    RealNumber hh_tag(analysis::UncertaintySource unc_source, analysis::UncertaintyScale unc_scale) const;
    Integer partonFlavour() const;
    Integer hadronFlavour() const;
    RealNumber rawf() const;
    RealNumber resolution() const;
    FilterBits triggerFilterMatch() const;
    size_t jet_index() const;

private:
    template<typename Value>
    Value CheckAndGet(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGet(jet_id, col, "jet", branch_name);
    }

    template<typename Value>
    const Value& CheckAndGetRef(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGetRef(jet_id, col, "jet", branch_name);
    }

private:
    size_t jet_id;
};

class TupleSubJet : public TupleObject {
public:
    TupleSubJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;

private:
    size_t jet_id;
};

class TupleFatJet : public TupleObject {
public:
    enum class MassType { Pruned, Filtered, Trimmed, SoftDrop };

    TupleFatJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    float m(MassType massType) const;
    DiscriminatorResult jettiness(size_t tau_index) const;
    const std::vector<TupleSubJet>& subJets() const;

private:
    template<typename Value>
    Value CheckAndGet(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGet(jet_id, col, "fatJet", branch_name);
    }

    template<typename Value>
    const Value& CheckAndGetRef(const std::vector<Value>& col, std::string_view branch_name) const
    {
        return TupleObject::CheckAndGetRef(jet_id, col, "fatJet", branch_name);
    }

private:
    size_t jet_id;
    std::vector<TupleSubJet> sub_jets;
};

class TupleMet : public TupleObject {
public:
    using CovMatrix = analysis::SquareMatrix<2>;
    TupleMet(const ntuple::Event& _event, MetType _met_type);
    MetType type() const;
    const LorentzVectorM& p4() const;
    const CovMatrix& cov() const;
    RealNumber pt() const;
    RealNumber phi() const;

private:
    MetType met_type;
};

} // namespace ntuple
