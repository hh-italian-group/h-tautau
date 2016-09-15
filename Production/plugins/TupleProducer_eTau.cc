/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::ETau;

    SelectionResults selection;
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) cut(triggerTools.HaveTriggerFired(hltPaths), "trigger");

    //Di-Lepton Veto
    const auto z_electrons = CollectZelectrons();
    const auto z_electrons_candidates = FindCompatibleObjects(z_electrons, z_electrons, DeltaR_DileptonVeto, "Z_e_e",0);
    selection.Zveto = z_electrons_candidates.size();

    // Signal-like leptons selection
    const auto selectedElectrons = CollectSignalElectrons();
    cut(selectedElectrons.size(),"electrons");
    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(),"taus");

    auto higgses = FindCompatibleObjects(selectedElectrons, selectedTaus, DeltaR_betweenSignalObjects, "H_e_tau");
    cut(higgses.size(),"ele_tau_pair");

    std::vector<HiggsCandidate> selected_higgses = higgses;
    if(applyTriggerMatch) {
        selected_higgses = triggerTools.ApplyTriggerMatch(higgses, hltPaths, false);
        cut(selected_higgses.size(), "triggerMatch");
    }

    std::sort(selected_higgses.begin(), selected_higgses.end(), &HiggsComparitor<HiggsCandidate>);
    selection.SetHiggsCandidate(selected_higgses.front());

    //Third-Lepton Veto
    const auto electronVetoCollection = CollectVetoElectrons(&selection.higgs->GetFirstDaughter());
    const auto muonVetoCollection = CollectVetoMuons();
    selection.electronVeto = electronVetoCollection.size();
    selection.muonVeto = muonVetoCollection.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer.Fit(*selection.higgs, *met);
    FillEventTuple(selection);
    FillSyncTuple(selection);
}

std::vector<BaseTupleProducer::ElectronCandidate> TupleProducer_eTau::CollectSignalElectrons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_eTau::SelectSignalElectron, this, _1, _2);
    return CollectObjects("SignalElectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_eTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_eTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_eTau::SelectSignalElectron(const ElectronCandidate& electron, Cutter& cut) const
{
    using namespace cuts::Htautau_2015::ETau;
    using namespace cuts::Htautau_2015::ETau::electronID;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electronD0 = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electronD0 < electronID::d0, "dxy", electronD0);
    const double electronDZ = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electronDZ < electronID::dz, "dz", electronDZ);
    const bool isTight = (*tight_id_decisions)[electron.getPtr()];
    cut(isTight, "electronMVATightID");
    const int eleMissingHits = electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    cut(eleMissingHits <= electronID::missingHits, "missingHits", eleMissingHits);
    cut(electron->passConversionVeto(), "conversionVeto");
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::Htautau_2015::ETau;
    using namespace cuts::Htautau_2015::ETau::tauID;

    cut(true, "gt0_tau_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    cut(p4.Pt() > pt, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const int dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > tauID::decayModeFinding, "oldDecayMode", dmFinding);
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == 1, "charge", tau->charge());
}

void TupleProducer_eTau::FillSyncTuple(const SelectionResults& selection)
{

  using namespace analysis;
  static const float default_value = ntuple::DefaultFillValue<Float_t>();

  BaseTupleProducer::FillSyncTuple(selection);
  syncTuple().pairType = static_cast<int>(analysis::Channel::ETau);

  // Leg 1, lepton
  const ElectronCandidate& electron = selection.higgs->GetFirstDaughter();
  syncTuple().dau1_pt  = electron.GetMomentum().Pt();
  syncTuple().dau1_eta  = electron.GetMomentum().Eta();
  syncTuple().dau1_phi  = electron.GetMomentum().Phi();
  syncTuple().dau1_iso   = electron.GetIsolation();

  syncTuple.Fill();
}

void TupleProducer_eTau::FillEventTuple(const SelectionResults& selection)
{
    using namespace analysis;
    static const float default_value = ntuple::DefaultFillValue<Float_t>();

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelID = static_cast<int>(analysis::Channel::ETau);

    // Leg 1, lepton
    const ElectronCandidate& electron = selection.higgs->GetFirstDaughter();
    eventTuple().p4_1     = analysis::LorentzVectorM(electron.GetMomentum());
    eventTuple().q_1      = electron.GetCharge();
    eventTuple().pfmt_1   = Calculate_MT(electron.GetMomentum(), met->GetMomentum().Pt(), met->GetMomentum().Phi());
    eventTuple().d0_1     = electron->gsfTrack()->dxy(primaryVertex->position());
    eventTuple().dZ_1     = electron->gsfTrack()->dz(primaryVertex->position());
    eventTuple().iso_1    = electron.GetIsolation();
    eventTuple().id_e_mva_nt_loose_1 = default_value;
    eventTuple().gen_match_1 = isMC ? gen_truth::genMatch(electron->p4(), *genParticles) : default_value;

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
