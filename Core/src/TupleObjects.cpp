/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include "h-tautau/Core/include/TupleObjects.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"

namespace ntuple {

TupleObject::TupleObject(const Event& _event) : event(&_event) {}

void TupleObject::CheckIndexRange(size_t index, size_t size, std::string_view obj_name,
                                  std::string_view branch_name) const
{
    if(index >= size) {
        analysis::EventIdentifier event_id(*event);
        throw analysis::exception("%1%: %2% index = %3% is out of range to index %4%, which has the size = %5%.")
                % event_id % obj_name % index % branch_name % size;
    }
}

TupleLepton::TupleLepton(const ntuple::Event& _event, size_t _object_id)
    : TupleObject(_event), object_id(_object_id)
{
}

const LorentzVectorM& TupleLepton::p4() const { return CheckAndGetRef(event->lep_p4, "lep_p4"); }
TupleObject::Integer TupleLepton::charge() const { return CheckAndGet(event->lep_q, "lep_q"); }
TupleObject::RealNumber TupleLepton::dxy() const { return CheckAndGet(event->lep_dxy, "lep_dxy"); }
TupleObject::RealNumber TupleLepton::dz() const { return CheckAndGet(event->lep_dz, "lep_dz"); }
TupleObject::RealNumber TupleLepton::iso() const { return CheckAndGet(event->lep_iso, "lep_iso"); }
analysis::GenLeptonMatch TupleLepton::gen_match() const
{
    return analysis::GenLeptonMatch(CheckAndGet(event->lep_gen_match, "lep_gen_match"));
}
const LorentzVectorM& TupleLepton::gen_p4() const { return CheckAndGetRef(event->lep_gen_p4, "lep_gen_p4"); }
TupleObject::Integer TupleLepton::decayMode() const { return CheckAndGet(event->lep_decayMode, "lep_decayMode"); }
analysis::LegType TupleLepton::leg_type() const
{
    return analysis::LegType(CheckAndGet(event->lep_type, "lep_type"));
}
bool TupleLepton::passConversionVeto() const
{
    return CheckAndGet(event->lep_elePassConversionVeto, "lep_elePassConversionVeto");
}
bool TupleLepton::passEleIsoId(DiscriminatorWP wp) const
{
    DiscriminatorIdResults eleId(CheckAndGet(event->lep_eleId_iso, "lep_eleId_iso"));
    return eleId.Passed(wp);
}
bool TupleLepton::passEleNoIsoId(DiscriminatorWP wp) const
{
    DiscriminatorIdResults eleId(CheckAndGet(event->lep_eleId_noIso, "lep_eleId_noIso"));
    return eleId.Passed(wp);
}
bool TupleLepton::passMuonId(DiscriminatorWP wp) const
{
    DiscriminatorIdResults muonId(CheckAndGet(event->lep_muonId, "lep_muonId"));
    return muonId.Passed(wp);
}

bool TupleLepton::Passed(analysis::TauIdDiscriminator tauIdDiscriminator, DiscriminatorWP wp) const
{
    if(leg_type() != analysis::LegType::tau)
        throw analysis::exception("LegType is not a tau in Passed for %1%") % analysis::EventIdentifier(*event);
    uint16_t discriminator_value = std::numeric_limits<uint16_t>::max();
    #define TAU_ID(name, pattern, has_raw, wp_list) \
        if(tauIdDiscriminator == analysis::TauIdDiscriminator::name) \
            discriminator_value = CheckAndGet(event->name, #name);
    TAU_IDS()
    #undef TAU_ID
    if(discriminator_value == std::numeric_limits<uint16_t>::max())
        throw analysis::exception("TauId discriminator value not found.");
    DiscriminatorIdResults discriminator(discriminator_value);
    return discriminator.Passed(wp);

}

bool TupleLepton::PassedOldDecayMode() const
{
    return CheckAndGet(event->lep_oldDecayModeFinding, "lep_oldDecayModeFinding");
}
bool TupleLepton::PassedNewDecayMode() const
{
    return CheckAndGet(event->lep_newDecayModeFinding, "lep_newDecayModeFinding");
}

TupleObject::DiscriminatorResult TupleLepton::GetRawValue(analysis::TauIdDiscriminator tauIdDiscriminator) const
{
    if(leg_type() != analysis::LegType::tau)
        throw analysis::exception("LegType is not a tau in Get Raw for %1%") % analysis::EventIdentifier(*event);
    #define TAU_ID(name, pattern, has_raw, wp_list) \
        if(tauIdDiscriminator == analysis::TauIdDiscriminator::name) \
            return CheckAndGet(event->name##raw, #name"raw");
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

TupleJet::TupleJet(const ntuple::Event& _event, size_t _jet_id) : TupleObject(_event), jet_id(_jet_id) {}
const LorentzVectorE& TupleJet::p4() const { return CheckAndGetRef(event->jets_p4, "jet_p4"); }
analysis::DiscriminatorIdResults TupleJet::GetPuId() const
{
    return analysis::DiscriminatorIdResults(CheckAndGet(event->jets_pu_id_upd, "jets_pu_id_upd"));
}

Float_t TupleJet::GetPuIdRaw() const { return CheckAndGet(event->jets_pu_id_upd_raw, "jets_pu_id_upd_raw"); }

bool TupleJet::PassPuId(DiscriminatorWP wp) const {
    const analysis::DiscriminatorIdResults jet_pu_id = GetPuId();
    return jet_pu_id.Passed(wp);
}
TupleObject::DiscriminatorResult TupleJet::csv() const { return CheckAndGet(event->jets_csv, "jets_csv"); }
TupleObject::DiscriminatorResult TupleJet::deepcsv() const
{
    return CheckAndGet(event->jets_deepCsv_BvsAll, "jets_deepCsv_BvsAll");
}
TupleObject::DiscriminatorResult TupleJet::deepFlavour() const
{
    return CheckAndGet(event->jets_deepFlavour_b, "jets_deepFlavour_b")
         + CheckAndGet(event->jets_deepFlavour_bb, "jets_deepFlavour_bb")
         + CheckAndGet(event->jets_deepFlavour_lepb, "jets_deepFlavour_lepb");
}
TupleObject::DiscriminatorResult TupleJet::deepFlavour_CvsL() const
{
    const auto prob_c = CheckAndGet(event->jets_deepFlavour_c, "jets_deepFlavour_c");
    const auto prob_uds = CheckAndGet(event->jets_deepFlavour_uds, "jets_deepFlavour_uds");
    const auto prob_g = CheckAndGet(event->jets_deepFlavour_g, "jets_deepFlavour_g");
    return prob_c / (prob_c + prob_uds + prob_g);
}
TupleObject::DiscriminatorResult TupleJet::deepFlavour_CvsB() const
{
    const auto prob_c = CheckAndGet(event->jets_deepFlavour_c, "jets_deepFlavour_c");
    const auto prob_b = CheckAndGet(event->jets_deepFlavour_b, "jets_deepFlavour_b");
    const auto prob_bb = CheckAndGet(event->jets_deepFlavour_bb, "jets_deepFlavour_bb");
    const auto prob_lepb = CheckAndGet(event->jets_deepFlavour_lepb, "jets_deepFlavour_lepb");
    return prob_c / (prob_c + prob_b + prob_bb + prob_lepb);
}
TupleObject::Integer TupleJet::partonFlavour() const
{
    return CheckAndGet(event->jets_partonFlavour, "jets_partonFlavour");
}
TupleObject::Integer TupleJet::hadronFlavour() const
{
    return CheckAndGet(event->jets_hadronFlavour, "jets_hadronFlavour");
}
TupleObject::RealNumber TupleJet::rawf() const { return CheckAndGet(event->jets_rawf, "jets_rawf"); }
TupleObject::RealNumber TupleJet::resolution() const { return CheckAndGet(event->jets_resolution, "jets_resolution"); }
size_t TupleJet::jet_index() const { return jet_id; }

TupleJet::FilterBits TupleJet::triggerFilterMatch() const
{
    analysis::TriggerDescriptorCollection::RootBitsContainer match_bits = {{
        CheckAndGet(event->jets_triggerFilterMatch_0, "jets_triggerFilterMatch_0"),
        CheckAndGet(event->jets_triggerFilterMatch_1, "jets_triggerFilterMatch_1"),
        CheckAndGet(event->jets_triggerFilterMatch_2, "jets_triggerFilterMatch_2"),
        CheckAndGet(event->jets_triggerFilterMatch_3, "jets_triggerFilterMatch_3"),
    }};
    return analysis::TriggerDescriptorCollection::ConvertFromRootRepresentation(match_bits);
}

TupleObject::RealNumber TupleJet::hh_btag() const
{
    if(!hh_btag_score.is_initialized())
        throw analysis::exception("TupleJet: hh_btag score is not set.");
    return *hh_btag_score;
}
void TupleJet::set_hh_btag(RealNumber score) { hh_btag_score = score; }

TupleSubJet::TupleSubJet(const ntuple::Event& _event, size_t _jet_id) : TupleObject(_event), jet_id(_jet_id) {}
const LorentzVectorE& TupleSubJet::p4() const
{
    return CheckAndGetRef(jet_id, event->subJets_p4, "sub-jet", "subJets_p4");
}

TupleFatJet::TupleFatJet(const ntuple::Event& _event, size_t _jet_id)
    : TupleObject(_event), jet_id(_jet_id)
{
    for(size_t n = 0; n < event->subJets_parentIndex.size(); ++n) {
        if(event->subJets_parentIndex.at(n) == jet_id)
            sub_jets.emplace_back(_event, n);
    }
}

const LorentzVectorE& TupleFatJet::p4() const { return CheckAndGetRef(event->fatJets_p4, "fatJets_p4"); }

float TupleFatJet::m(MassType massType) const
{
    if(massType == MassType::SoftDrop) return CheckAndGet(event->fatJets_m_softDrop, "fatJets_m_softDrop");
    throw analysis::exception("Unsupported fat jet mass type");
}

TupleObject::DiscriminatorResult TupleFatJet::jettiness(size_t tau_index) const
{
    if(tau_index == 1) return CheckAndGet(event->fatJets_jettiness_tau1, "fatJets_jettiness_tau1");
    if(tau_index == 2) return CheckAndGet(event->fatJets_jettiness_tau2, "fatJets_jettiness_tau2");
    if(tau_index == 3) return CheckAndGet(event->fatJets_jettiness_tau3, "fatJets_jettiness_tau3");
    if(tau_index == 4) return CheckAndGet(event->fatJets_jettiness_tau4, "fatJets_jettiness_tau4");
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
const LorentzVectorM& TupleMet::p4() const { return event->pfMET_p4; }
const TupleMet::CovMatrix& TupleMet::cov() const { return event->pfMET_cov; }
TupleObject::RealNumber TupleMet::pt() const { return p4().pt(); }
TupleObject::RealNumber TupleMet::phi() const { return p4().phi(); }

} // namespace ntuple
