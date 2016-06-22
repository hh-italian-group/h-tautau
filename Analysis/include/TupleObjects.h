/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using LorentzVectorE = analysis::LorentzVector;
    using LorentzVectorM = analysis::LorentzVectorM;
    using DiscriminatorResult = float;

    TupleObject(const ntuple::Event& _event) : event(_event) {}

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

        momentum = leg_id == 1
                 ? LorentzVectorM(event->pt_1, event->eta_1, event->phi_1, event->m_1)
                 : LorentzVectorM(event->pt_2, event->eta_2, event->phi_2, event->m_2);
    }

    const LorentzVectorM& p4() const { return momentum; }
    int charge() const { return leg_id == 1 ? event->q_1 : event->q_2; }
    double d0() const { return leg_id == 1 ? event->d0_1 : event->d0_2; }
    double dZ() const { return leg_id == 1 ? event->dZ_1 : event->dZ_2; }

    double mt(MetType metType = MetType::PF) const
    {
        if(metType == MetType::PF)
            return leg_id == 1 ? event->pfmt_1 : event->pfmt_2;
        if(metType == MetType::MVA)
            return leg_id == 1 ? event->mt_1 : event->mt_2;
        if(metType == MetType::PUPPI)
            return leg_id == 1 ? event->puppimt_1 : event->puppimt_2;
        throw analysis::exception("Unsupported MET type.");
    }

    double iso() const { return leg_id == 1 ? event->iso_1 : event->iso_2; }
    int gen_match() const { return leg_id == 1 ? event->gen_match_1 : event->gen_match_2; }

protected:
    size_t leg_id;
    LorentzVectorM momentum;
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
    using TupleLepton::TupleLepton;

    DiscriminatorResult againstElectronMVA6(DiscriminatorWP wp) const
    {
        if(wp == DiscriminatorWP::VLoose)
            return leg_id == 1 ? event->againstElectronVLooseMVA6_1 : event->againstElectronVLooseMVA6_2;
        if(wp == DiscriminatorWP::Loose)
            return leg_id == 1 ? event->againstElectronLooseMVA6_1 : event->againstElectronLooseMVA6_2;
        if(wp == DiscriminatorWP::Medium)
            return leg_id == 1 ? event->againstElectronMediumMVA6_1 : event->againstElectronMediumMVA6_2;
        if(wp == DiscriminatorWP::Tight)
            return leg_id == 1 ? event->againstElectronTightMVA6_1 : event->againstElectronTightMVA6_2;
        if(wp == DiscriminatorWP::VTight)
            return leg_id == 1 ? event->againstElectronVTightMVA6_1 : event->againstElectronVTightMVA6_2;
        throw analysis::exception("againstElectronMVA6: unsupported discriminator WP.");
    }

    DiscriminatorResult byIsolationMVA2raw() const
    {
        return leg_id == 1 ? event->byIsolationMVA2raw_1 : event->byIsolationMVA2raw_2;
    }

    DiscriminatorResult againstMuon3(DiscriminatorWP wp) const
    {
        if(wp == DiscriminatorWP::Loose)
            return leg_id == 1 ? event->againstMuonLoose3_1 : event->againstMuonLoose3_2;
        if(wp == DiscriminatorWP::Tight)
            return leg_id == 1 ? event->againstMuonTight3_1 : event->againstMuonTight3_2;
        throw analysis::exception("againstMuon3: unsupported discriminator WP.");
    }

    DiscriminatorResult byCombinedIsolationDeltaBetaCorrRaw3Hits() const
    {
        return leg_id == 1 ? event->byCombinedIsolationDeltaBetaCorrRaw3Hits_1
                           : event->byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
    }

    DiscriminatorResult byIsolationMVA3raw(bool use_new_dm, bool use_lifetime) const
    {
        if(use_new_dm) {
            if(use_lifetime)
                return leg_id == 1 ? event->byIsolationMVA3newDMwLTraw_1 : event->byIsolationMVA3newDMwLTraw_2;
            return leg_id == 1 ? event->byIsolationMVA3newDMwoLTraw_1 : event->byIsolationMVA3newDMwoLTraw_2;
        }
        if(use_lifetime)
            return leg_id == 1 ? event->byIsolationMVA3oldDMwLTraw_1 : event->byIsolationMVA3oldDMwLTraw_2;
        return leg_id == 1 ? event->byIsolationMVA3oldDMwoLTraw_1 : event->byIsolationMVA3oldDMwoLTraw_2;
    }
};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(jet_id)
    {
        if(jet_id >= event->pt_jets.size())
            throw analysis::exception("Jet id = %1% is out of range.") % jet_id;

        momentum = LorentzVectorE(event->pt_jets.at(jet_id), event->eta_jets.at(jet_id),
                                  event->phi_jets.at(jet_id), event->energy_jets.at(jet_id));
    }

    const LorentzVectorE& p4() const { return momentum; }
    DiscriminatorResult mva() const { return event->mva_jets.at(jet_id); }
    DiscriminatorResult csv() const { return event->csv_jets.at(jet_id); }
    DiscriminatorResult partonFlavour() const { return event->partonFlavour_jets.at(jet_id); }

protected:
    size_t jet_id;
    LorentzVectorE momentum;
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

        if(met_type == MetType::PF) {
            cov_matrix[0][0] = event->metcov00;
            cov_matrix[0][1] = event->metcov01;
            cov_matrix[1][0] = event->metcov10;
            cov_matrix[1][1] = event->metcov11;
        }
        if(met_type == MetType::MVA) {
            cov_matrix[0][0] = event->mvacov00;
            cov_matrix[0][1] = event->mvacov01;
            cov_matrix[1][0] = event->mvacov10;
            cov_matrix[1][1] = event->mvacov11;
        }
    }

    MetType type() const { return met_type; }

    double pt() const
    {
        if(met_type == MetType::PF) return event->met;
        if(met_type == MetType::MVA) return event->mvamet;
        return event->puppimet;
    }

    double phi() const
    {
        if(met_type == MetType::PF) return event->metphi;
        if(met_type == MetType::MVA) return event->mvametphi;
        return event->puppimetphi;
    }

    const CovMatrix& cov() const { return cov_matrix; }

protected:
    MetType met_type;
    CovMatrix cov_matrix;
};

} // namespace ntuple
