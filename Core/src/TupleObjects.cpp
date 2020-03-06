/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include "h-tautau/Core/include/TupleObjects.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"

namespace ntuple {

TupleObject::TupleObject(const ntuple::Event& _event) : event(&_event) {}

TupleLepton::TupleLepton(const ntuple::Event& _event, size_t _object_id)
    : TupleObject(_event), object_id(_object_id)
{
    if(object_id >= event->lep_p4.size())
        throw analysis::exception("Invalid object id = %1%.") % object_id;
}

const LorentzVectorM& TupleLepton::p4() const { return event->lep_p4.at(object_id); }
TupleObject::Integer TupleLepton::charge() const { return event->lep_q.at(object_id); }
TupleObject::RealNumber TupleLepton::dxy() const { return event->lep_dxy.at(object_id); }
TupleObject::RealNumber TupleLepton::dz() const { return event->lep_dz.at(object_id); }
TupleObject::RealNumber TupleLepton::iso() const { return event->lep_iso.at(object_id); }
analysis::GenLeptonMatch TupleLepton::gen_match() const { return analysis::GenLeptonMatch(event->lep_gen_match.at(object_id)); }
const LorentzVectorM& TupleLepton::gen_p4() const { return event->lep_gen_p4.at(object_id); }
TupleObject::Integer TupleLepton::decayMode() const { return event->lep_decayMode.at(object_id); }
analysis::LegType TupleLepton::leg_type() const { return analysis::LegType(event->lep_type.at(object_id)); }
bool TupleLepton::passConversionVeto() const {return event->lep_elePassConversionVeto.at(object_id);}
bool TupleLepton::passEleIso(DiscriminatorWP wp) const
{
    DiscriminatorIdResults eleIso(event->lep_eleId_iso.at(object_id));
    return eleIso.Passed(wp);
}

bool TupleLepton::passMuonId(DiscriminatorWP wp) const
{
    DiscriminatorIdResults muonId(event->lep_muonId.at(object_id));
    return muonId.Passed(wp);
}

bool TupleLepton::Passed(analysis::TauIdDiscriminator tauIdDiscriminator, DiscriminatorWP wp) const
{
    if(leg_type() != analysis::LegType::tau)
        throw analysis::exception("LegType is not a tau in Passed for %1%") % analysis::EventIdentifier(*event);
    uint16_t discriminator_value = std::numeric_limits<uint16_t>::max();
    #define TAU_ID(name, pattern, has_raw, wp_list) \
        if(tauIdDiscriminator == analysis::TauIdDiscriminator::name) discriminator_value = event->name.at(object_id);
    TAU_IDS()
    #undef TAU_ID
    if(discriminator_value == std::numeric_limits<uint16_t>::max())
        throw analysis::exception("TauId discriminator value not found.");
    DiscriminatorIdResults discriminator(discriminator_value);
    return discriminator.Passed(wp);

}

bool TupleLepton::PassedOldDecayMode() const { return event->lep_oldDecayModeFinding.at(object_id); }
bool TupleLepton::PassedNewDecayMode() const { return event->lep_newDecayModeFinding.at(object_id); }

TupleObject::DiscriminatorResult TupleLepton::GetRawValue(analysis::TauIdDiscriminator tauIdDiscriminator) const
{
    if(leg_type() != analysis::LegType::tau)
        throw analysis::exception("LegType is not a tau in Get Raw for %1%") % analysis::EventIdentifier(*event);
    #define TAU_ID(name, pattern, has_raw, wp_list) \
        if(tauIdDiscriminator == analysis::TauIdDiscriminator::name) return event->name##raw.at(object_id);
    TAU_IDS()
    #undef TAU_ID
    throw analysis::exception("TauId Raw value not found.");
}

int TupleLepton::CompareIsolations(const TupleLepton& other, analysis::TauIdDiscriminator disc) const
{
    if(leg_type() != other.leg_type())
        throw analysis::exception("Isolation of legs with different types are not comparable");
    if(leg_type() == analysis::LegType::e || leg_type() == analysis::LegType::mu) {
        if(iso() == other.iso()) return 0;
        return iso() < other.iso() ? 1 : -1;
    }
    if(leg_type() == analysis::LegType::tau) {
        const auto iso1 = GetRawValue(disc);
        const auto iso2 = other.GetRawValue(disc);
        if(iso1 == iso2) return 0;
        return iso1 > iso2 ? 1 : -1;
    }
    throw analysis::exception("Isolation comparison for the leg type '%1%' is not supported.") % leg_type();
}

TupleJet::TupleJet(const ntuple::Event& _event, size_t _jet_id)
    : TupleObject(_event), jet_id(_jet_id)
{
    if(jet_id >= event->jets_p4.size())
        throw analysis::exception("Jet id = %1% is out of range.") % jet_id;
}

const LorentzVectorE& TupleJet::p4() const { return event->jets_p4.at(jet_id); }
analysis::DiscriminatorIdResults TupleJet::GetPuId() const
{
    analysis::DiscriminatorIdResults jet_pu_id(event->jets_pu_id_upd.at(jet_id));
    return jet_pu_id;
}

Float_t TupleJet::GetPuIdRaw() const
{
    return event->jets_pu_id_upd_raw.at(jet_id);
}

bool TupleJet::PassPuId(DiscriminatorWP wp) const {
    analysis::DiscriminatorIdResults jet_pu_id = GetPuId();
    return jet_pu_id.Passed(wp);
}
TupleObject::DiscriminatorResult TupleJet::csv() const { return event->jets_csv.at(jet_id); }
TupleObject::DiscriminatorResult TupleJet::deepcsv() const { return event->jets_deepCsv_BvsAll.at(jet_id); }
TupleObject::DiscriminatorResult TupleJet::deepFlavour() const
{
    return event->jets_deepFlavour_b.at(jet_id) + event->jets_deepFlavour_bb.at(jet_id)
           + event->jets_deepFlavour_lepb.at(jet_id);
}
TupleObject::RealNumber TupleJet::hh_tag(analysis::UncertaintySource unc_source,
                                         analysis::UncertaintyScale unc_scale) const
{
    for(unsigned n = 0; n < event->jet_hh_score_value.size(); ++n){
        if(event->jet_hh_score_index.at(n) != jet_id) continue;
        if(event->jet_hh_score_unc_source.at(n) != static_cast<int>(unc_source)) continue;
        if(event->jet_hh_score_unc_scale.at(n) != static_cast<int>(unc_scale)) continue;
        return event->jet_hh_score_value.at(n);
    }
    throw analysis::exception("HH tag not found for TupleObjects.");
}
TupleObject::Integer TupleJet::partonFlavour() const { return event->jets_partonFlavour.at(jet_id); }
TupleObject::Integer TupleJet::hadronFlavour() const { return event->jets_hadronFlavour.at(jet_id); }
TupleObject::RealNumber TupleJet::rawf() const { return event->jets_rawf.at(jet_id); }
TupleObject::RealNumber TupleJet::resolution() const { return event->jets_resolution.at(jet_id); }
size_t TupleJet::jet_index() const { return jet_id; }

TupleJet::FilterBits TupleJet::triggerFilterMatch() const
{
    analysis::TriggerDescriptorCollection::RootBitsContainer match_bits = {{
        event->jets_triggerFilterMatch_0.at(jet_id),
        event->jets_triggerFilterMatch_1.at(jet_id),
        event->jets_triggerFilterMatch_2.at(jet_id),
        event->jets_triggerFilterMatch_3.at(jet_id),
    }};
    return analysis::TriggerDescriptorCollection::ConvertFromRootRepresentation(match_bits);
}

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
