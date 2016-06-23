/*! Implementation of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_tauTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_tauTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::MuTau;

    SelectionResults selection;
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) cut(triggerTools.HaveTriggerFired(hltPaths), "trigger");

    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");

    auto higgses = FindCompatibleObjects(selectedTaus, selectedTaus, DeltaR_betweenSignalObjects, "H_tau_tau");
    cut(higgses.size(), "tau_tau_pair");

    std::vector<HiggsCandidate> selected_higgses = higgses;
    if(applyTriggerMatch) {
        auto triggeredHiggses = triggerTools.ApplyTriggerMatch(higgses, hltPaths, true);
        cut(triggeredHiggses.size(), "triggerMatch");

        selected_higgses = triggerTools.ApplyL1TriggerTauMatch(triggeredHiggses);
        cut(selected_higgses.size(), "L1triggerMatch");
    }

    std::sort(selected_higgses.begin(), selected_higgses.end(), &HiggsComparitor<HiggsCandidate>);
    selection.SetHiggsCandidate(selected_higgses.front());

    //Third-Lepton Veto
    const auto electronVetoCollection = CollectVetoElectrons();
    const auto muonVetoCollection = CollectVetoMuons();
    selection.electronVeto = electronVetoCollection.size();
    selection.muonVeto = muonVetoCollection.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer.Fit(*selection.higgs, *met);
    FillEventTuple(selection);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_tauTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_tauTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_tauTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::Htautau_2015::TauTau;
    using namespace cuts::Htautau_2015::TauTau::tauID;

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

void TupleProducer_tauTau::FillEventTuple(const SelectionResults& selection)
{
    using namespace analysis;
    static const float default_value = ntuple::DefaultFillValue<Float_t>();

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelID = static_cast<int>(analysis::Channel::TauTau);

    // Leg 1, tau
    const TauCandidate& tau = selection.higgs->GetFirstDaughter();
    eventTuple().p4_1     = analysis::LorentzVectorM(tau.GetMomentum());
    eventTuple().q_1      = tau.GetCharge();
    eventTuple().pfmt_1   = Calculate_MT(tau.GetMomentum(), met->GetMomentum().Pt(), met->GetMomentum().Phi());
    eventTuple().d0_1     = Calculate_dxy(tau->vertex(), primaryVertex->position(), tau.GetMomentum());
    eventTuple().dZ_1     = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())->dz();
    eventTuple().iso_1    = tau.GetIsolation();
    eventTuple().gen_match_1 = isMC ? gen_truth::genMatch(tau->p4(), *genParticles) : default_value;
    eventTuple().tauIDs_1.insert(tau->tauIDs().begin(), tau->tauIDs().end());

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau);
