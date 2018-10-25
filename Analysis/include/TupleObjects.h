/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"
#include "TauIdResults.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using DiscriminatorResult = float;
    using Integer = int;
    using RealNumber = float;

    TupleObject(const ntuple::Event& _event) : event(&_event) {}

protected:
    const Event* event;
};

class TupleLepton : public TupleObject {
public:
    TupleLepton(const ntuple::Event& _event, size_t _leg_id)
        : TupleObject(_event), leg_id(_leg_id)
    {
        if(leg_id < 1 || leg_id > 2)
            throw analysis::exception("Invalid leg id = %1%.") % leg_id;
    }

    const LorentzVectorM& p4() const { return leg_id == 1 ? event->p4_1 : event->p4_2; }
    Integer charge() const { return leg_id == 1 ? event->q_1 : event->q_2; }
    RealNumber dxy() const { return leg_id == 1 ? event->dxy_1 : event->dxy_2; }
    RealNumber dz() const { return leg_id == 1 ? event->dz_1 : event->dz_2; }
    RealNumber iso() const { return leg_id == 1 ? event->iso_1 : event->iso_2; }
    Integer gen_match() const { return leg_id == 1 ? event->gen_match_1 : event->gen_match_2; }

protected:
    size_t leg_id;
};

class TupleElectron : public TupleLepton {
public:
    explicit TupleElectron(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleMuon : public TupleLepton {
public:
    explicit TupleMuon(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleTau : public TupleLepton {
public:
    TupleTau(const ntuple::Event& _event, size_t _leg_id)
        : TupleLepton(_event, _leg_id), tauIds(leg_id == 1 ? event->tauId_flags_1 : event->tauId_flags_2)
    {
    }

public:
    const analysis::TauIdResults& tauIDs() const { return tauIds; }

    bool tauID(analysis::TauIdDiscriminator discriminator, analysis::DiscriminatorWP wp) const
    {
        return tauIds.Result(discriminator, wp);
    }

    DiscriminatorResult tauIDraw(analysis::TauIdDiscriminator discriminator) const
    {
        using RawValuePtr = const DiscriminatorResult Event::*;
        using TauIdDiscriminator = analysis::TauIdDiscriminator;
        using Key = std::pair<TauIdDiscriminator, size_t>;
        static const std::map<Key, RawValuePtr> raw_values = {
            { { TauIdDiscriminator::againstElectronMVA6, 1 }, &Event::tauId_againstElectronMVA6Raw_1 },
            { { TauIdDiscriminator::againstElectronMVA6, 2 }, &Event::tauId_againstElectronMVA6Raw_2 },
            { { TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits, 1 },
                &Event::tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_1 },
            { { TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits, 2 },
                &Event::tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_2 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT, 1 },
                &Event::tauId_byIsolationMVArun2v1DBoldDMwLTraw_1 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT, 2 },
                &Event::tauId_byIsolationMVArun2v1DBoldDMwLTraw_2 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBdR03oldDMwLT, 1 },
                &Event::tauId_byIsolationMVArun2v1DBdR03oldDMwLTraw_1 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBdR03oldDMwLT, 2 },
                &Event::tauId_byIsolationMVArun2v1DBdR03oldDMwLTraw_2 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016, 1 },
                &Event::tauId_byIsolationMVArun2v1DBoldDMwLTraw2016_1 },
            { { TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016, 2 },
                &Event::tauId_byIsolationMVArun2v1DBoldDMwLTraw2016_2 },
            { { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, 1 },
                &Event::tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_1 },
            { { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, 2 },
                &Event::tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_2 },
            { { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, 1 },
                &Event::tauId_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017_1 },
            { { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, 2 },
                &Event::tauId_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017_2 }
        };
        const Key key(discriminator, leg_id);
        auto iter = raw_values.find(key);
        if(iter == raw_values.end())
            throw analysis::exception("Raw value not found for tau ID discriminator '%1%'.") % discriminator;
        return event->*(iter->second);
    }

private:
    analysis::TauIdResults tauIds;
};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->jets_p4.size())
            throw analysis::exception("Jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->jets_p4.at(jet_id); }
    bool PassPuId(DiscriminatorWP wp) const {
        if(wp == DiscriminatorWP::Loose)
            return event->jets_pu_id.at(jet_id) & (1 << 2);
        if(wp == DiscriminatorWP::Medium)
            return event->jets_pu_id.at(jet_id) & (1 << 1);
        if(wp == DiscriminatorWP::Tight)
            return event->jets_pu_id.at(jet_id) & (1 << 0);
        return false;
    }
    DiscriminatorResult csv() const { return event->jets_csv.at(jet_id); }
    DiscriminatorResult deepcsv() const { return event->jets_deepCsv_BvsAll.at(jet_id); }
    DiscriminatorResult deepFlavour() const{return (event->jets_deepFlavour_b.at(jet_id) +
                                event->jets_deepFlavour_bb.at(jet_id)+ event->jets_deepFlavour_lepb.at(jet_id))}
    Integer hadronFlavour() const { return event->jets_hadronFlavour.at(jet_id); }
    RealNumber rawf() const { return event->jets_rawf.at(jet_id); }
    RealNumber resolution() const { return event->jets_resolution.at(jet_id); }
    ULong64_t triggerFilterMatch() const { return event->jets_triggerFilterMatch.at(jet_id); } 
private:
    size_t jet_id;
};

class TupleSubJet : public TupleObject {
public:
    TupleSubJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->subJets_p4.size())
            throw analysis::exception("Fat sub-jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->subJets_p4.at(jet_id); }

private:
    size_t jet_id;
};

class TupleFatJet : public TupleObject {
public:
    enum class MassType { Pruned, Filtered, Trimmed, SoftDrop };

    TupleFatJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->fatJets_p4.size())
            throw analysis::exception("Fat jet id = %1% is out of range.") % jet_id;

        for(size_t n = 0; n < event->subJets_p4.size(); ++n) {
            if(event->subJets_parentIndex.at(n) == jet_id)
                sub_jets.push_back(TupleSubJet(_event, n));
        }
    }

    const LorentzVectorE& p4() const { return event->fatJets_p4.at(jet_id); }

    float m(MassType massType) const
    {
        if(massType == MassType::SoftDrop) return event->fatJets_m_softDrop.at(jet_id);
        throw analysis::exception("Unsupported fat jet mass type");
    }

    DiscriminatorResult jettiness(size_t tau_index) const
    {
        if(tau_index == 1) return event->fatJets_jettiness_tau1.at(jet_id);
        if(tau_index == 2) return event->fatJets_jettiness_tau2.at(jet_id);
        if(tau_index == 3) return event->fatJets_jettiness_tau3.at(jet_id);
        if(tau_index == 4) return event->fatJets_jettiness_tau4.at(jet_id);
        throw analysis::exception("Unsupported tau index = %1% for fat jet subjettiness.") % tau_index;
    }

    const std::vector<TupleSubJet>& subJets() const { return sub_jets; }

private:
    size_t jet_id;
    std::vector<TupleSubJet> sub_jets;
};

class TupleMet : public TupleObject {
public:
    using CovMatrix = analysis::SquareMatrix<2>;
    TupleMet(const ntuple::Event& _event, MetType _met_type)
        : TupleObject(_event), met_type(_met_type)
    {
        static const std::set<MetType> supported_types = { MetType::PF, MetType::MVA, MetType::PUPPI };
        if(!supported_types.count(met_type))
            throw analysis::exception("Unsupported met type.");
    }

    MetType type() const { return met_type; }

    const LorentzVectorM& p4() const
    {
        return event->pfMET_p4;
    }

    const CovMatrix& cov() const
    {
        return event->pfMET_cov;
    }

    RealNumber pt() const { return p4().pt(); }
    RealNumber phi() const { return p4().phi(); }

private:
    MetType met_type;
};

} // namespace ntuple
