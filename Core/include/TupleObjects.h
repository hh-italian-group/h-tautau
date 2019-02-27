/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"
#include "DiscriminatorIdResults.h"
#include "TauIdResults.h"

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
};

class TupleLepton : public TupleObject {
public:
    TupleLepton(const ntuple::Event& _event, size_t _object_id);
    const LorentzVectorM& p4() const;
    Integer charge() const;
    RealNumber dxy() const;
    RealNumber dz() const;
    RealNumber iso() const;
    Integer gen_match() const;
    Integer decayMode() const;

protected:
    size_t object_id;
};

class TupleElectron : public TupleLepton {
public:
    explicit TupleElectron(const ntuple::Event& _event, size_t _object_id = 0);
    DiscriminatorIdResults iso_wp() const { return DiscriminatorIdResults(); }
};

class TupleMuon : public TupleLepton {
public:
    explicit TupleMuon(const ntuple::Event& _event, size_t _object_id = 0);
    DiscriminatorIdResults iso_wp() const { return DiscriminatorIdResults(); }
};

class TupleTau : public TupleLepton {
public:
    using IdKey = uint32_t;
    using TupleLepton::TupleLepton;
    using ValueKeyPair = std::pair<std::string, IdKey>;

public:
    explicit TupleTau(const ntuple::Event& _event, size_t _object_id = 0);
    bool Passed(analysis::TauIdDiscriminator tauIdDiscriminator, DiscriminatorWP wp) const;
    Float_t GetRawValue(analysis::TauIdDiscriminator tauIdDiscriminator) const;

};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    bool PassPuId(DiscriminatorWP wp) const;
    DiscriminatorResult csv() const;
    DiscriminatorResult deepcsv() const;
    DiscriminatorResult deepFlavour() const;
    Integer hadronFlavour() const;
    RealNumber rawf() const;
    RealNumber resolution() const;
    ULong64_t triggerFilterMatch() const;
    size_t jet_index() const;
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
