/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include "h-tautau/Core/include/TupleObjects.h"

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TauIdResults.h"

namespace ntuple {

TupleObject::TupleObject(const ntuple::Event& _event) : event(&_event) {}

TupleLepton::TupleLepton(const ntuple::Event& _event, size_t _leg_id)
    : TupleObject(_event), leg_id(_leg_id)
{
    if(leg_id < 1 || leg_id > 2)
        throw analysis::exception("Invalid leg id = %1%.") % leg_id;
}

const LorentzVectorM& TupleLepton::p4() const { return leg_id == 1 ? event->p4_1 : event->p4_2; }
TupleObject::Integer TupleLepton::charge() const { return leg_id == 1 ? event->q_1 : event->q_2; }
TupleObject::RealNumber TupleLepton::dxy() const { return leg_id == 1 ? event->dxy_1 : event->dxy_2; }
TupleObject::RealNumber TupleLepton::dz() const { return leg_id == 1 ? event->dz_1 : event->dz_2; }
TupleObject::RealNumber TupleLepton::iso() const { return leg_id == 1 ? event->iso_1 : event->iso_2; }
TupleObject::Integer TupleLepton::gen_match() const { return leg_id == 1 ? event->gen_match_1 : event->gen_match_2; }

TupleElectron::TupleElectron(const ntuple::Event& _event, size_t _leg_id) : TupleLepton(_event, _leg_id) {}
TupleMuon::TupleMuon(const ntuple::Event& _event, size_t _leg_id) : TupleLepton(_event, _leg_id) {}

TupleTau::TupleTau(const ntuple::Event& _event, size_t _leg_id)
    : TupleLepton(_event, _leg_id), tauIds(leg_id == 1 ? event->tauId_flags_1 : event->tauId_flags_2)
{
}

const analysis::TauIdResults& TupleTau::tauIDs() const { return tauIds; }

bool TupleTau::tauID(analysis::TauIdDiscriminator discriminator, analysis::DiscriminatorWP wp) const
{
    return tauIds.Result(discriminator, wp);
}

TupleObject::DiscriminatorResult TupleTau::tauIDraw(analysis::TauIdDiscriminator discriminator) const
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

TupleJet::TupleJet(const ntuple::Event& _event, size_t _jet_id)
    : TupleObject(_event), jet_id(_jet_id)
{
    if(jet_id >= event->jets_p4.size())
        throw analysis::exception("Jet id = %1% is out of range.") % jet_id;
}

const LorentzVectorE& TupleJet::p4() const { return event->jets_p4.at(jet_id); }
bool TupleJet::PassPuId(DiscriminatorWP wp) const {
    if(wp == DiscriminatorWP::Loose)
        return event->jets_pu_id.at(jet_id) & (1 << 2);
    if(wp == DiscriminatorWP::Medium)
        return event->jets_pu_id.at(jet_id) & (1 << 1);
    if(wp == DiscriminatorWP::Tight)
        return event->jets_pu_id.at(jet_id) & (1 << 0);
    return false;
}
TupleObject::DiscriminatorResult TupleJet::csv() const { return event->jets_csv.at(jet_id); }
TupleObject::DiscriminatorResult TupleJet::deepcsv() const { return event->jets_deepCsv_BvsAll.at(jet_id); }
TupleObject::DiscriminatorResult TupleJet::deepFlavour() const
{
    return event->jets_deepFlavour_b.at(jet_id) + event->jets_deepFlavour_bb.at(jet_id)
           + event->jets_deepFlavour_lepb.at(jet_id);
}
TupleObject::Integer TupleJet::hadronFlavour() const { return event->jets_hadronFlavour.at(jet_id); }
TupleObject::RealNumber TupleJet::rawf() const { return event->jets_rawf.at(jet_id); }
TupleObject::RealNumber TupleJet::resolution() const { return event->jets_resolution.at(jet_id); }
ULong64_t TupleJet::triggerFilterMatch() const { return event->jets_triggerFilterMatch.at(jet_id); }
size_t TupleJet::jet_index() const { return jet_id; }

TupleSubJet::TupleSubJet(const ntuple::Event& _event, size_t _jet_id)
    : TupleObject(_event), jet_id(_jet_id)
{
    if(jet_id >= event->subJets_p4.size())
        throw analysis::exception("Fat sub-jet id = %1% is out of range.") % jet_id;
}

const LorentzVectorE& TupleSubJet::p4() const { return event->subJets_p4.at(jet_id); }

TupleFatJet::TupleFatJet(const ntuple::Event& _event, size_t _jet_id)
    : TupleObject(_event), jet_id(_jet_id)
{
    if(jet_id >= event->fatJets_p4.size())
        throw analysis::exception("Fat jet id = %1% is out of range.") % jet_id;

    for(size_t n = 0; n < event->subJets_p4.size(); ++n) {
        if(event->subJets_parentIndex.at(n) == jet_id)
            sub_jets.push_back(TupleSubJet(_event, n));
    }
}

const LorentzVectorE& TupleFatJet::p4() const { return event->fatJets_p4.at(jet_id); }

float TupleFatJet::m(MassType massType) const
{
    if(massType == MassType::SoftDrop) return event->fatJets_m_softDrop.at(jet_id);
    throw analysis::exception("Unsupported fat jet mass type");
}

TupleObject::DiscriminatorResult TupleFatJet::jettiness(size_t tau_index) const
{
    if(tau_index == 1) return event->fatJets_jettiness_tau1.at(jet_id);
    if(tau_index == 2) return event->fatJets_jettiness_tau2.at(jet_id);
    if(tau_index == 3) return event->fatJets_jettiness_tau3.at(jet_id);
    if(tau_index == 4) return event->fatJets_jettiness_tau4.at(jet_id);
    throw analysis::exception("Unsupported tau index = %1% for fat jet subjettiness.") % tau_index;
}

const std::vector<TupleSubJet>& TupleFatJet::subJets() const { return sub_jets; }

TupleMet::TupleMet(const ntuple::Event& _event, MetType _met_type)
    : TupleObject(_event), met_type(_met_type)
{
    static const std::set<MetType> supported_types = { MetType::PF, MetType::MVA, MetType::PUPPI };
    if(!supported_types.count(met_type))
        throw analysis::exception("Unsupported met type.");
}

TupleMet::MetType TupleMet::type() const { return met_type; }

const LorentzVectorM& TupleMet::p4() const
{
    return event->pfMET_p4;
}

const TupleMet::CovMatrix& TupleMet::cov() const
{
    return event->pfMET_cov;
}

TupleObject::RealNumber TupleMet::pt() const { return p4().pt(); }
TupleObject::RealNumber TupleMet::phi() const { return p4().phi(); }

} // namespace ntuple
