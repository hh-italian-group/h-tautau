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

    bool Passed(analysis::TauIdDiscriminator tauIdDiscriminator, DiscriminatorWP wp) const;
    bool PassedOldDecayMode() const;
    DiscriminatorResult GetRawValue(analysis::TauIdDiscriminator tauIdDiscriminator) const;

    int CompareIsolations(const TupleLepton& other, analysis::TauIdDiscriminator disc) const;

private:
    size_t object_id;
};


class TupleJet : public TupleObject {
public:
    using FilterBits = analysis::TriggerDescriptorCollection::BitsContainer;

    TupleJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    bool PassPuId(DiscriminatorWP wp) const;
    analysis::DiscriminatorIdResults GetPuId() const
    DiscriminatorResult csv() const;
    DiscriminatorResult deepcsv() const;
    DiscriminatorResult deepFlavour() const;
    Integer partonFlavour() const;
    Integer hadronFlavour() const;
    RealNumber rawf() const;
    RealNumber resolution() const;
    FilterBits triggerFilterMatch() const;
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
