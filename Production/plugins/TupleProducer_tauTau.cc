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

    cut(triggerTools.HaveTriggerFired(hltPaths), "trigger");

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

    eventTuple().pt_1     = tau.GetMomentum().Pt();
    eventTuple().phi_1    = tau.GetMomentum().Phi();
    eventTuple().eta_1    = tau.GetMomentum().Eta();
    eventTuple().m_1      = tau.GetMomentum().M();
    eventTuple().q_1      = tau.GetCharge();
    eventTuple().pfmt_1   = Calculate_MT(tau.GetMomentum(), met->GetMomentum().Pt(), met->GetMomentum().Phi());
    eventTuple().d0_1     = Calculate_dxy(tau->vertex(), primaryVertex->position(), tau.GetMomentum());
    eventTuple().dZ_1     = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())->dz();
    eventTuple().iso_1    = tau.GetIsolation();
    eventTuple().gen_match_2 = isMC ? gen_truth::genMatch(tau->p4(), *genParticles) : default_value;

    eventTuple().againstElectronLooseMVA6_1   = tau->tauID("againstElectronLooseMVA6");
    eventTuple().againstElectronMediumMVA6_1  = tau->tauID("againstElectronMediumMVA6");
    eventTuple().againstElectronTightMVA6_1   = tau->tauID("againstElectronTightMVA6");
    eventTuple().againstElectronVLooseMVA6_1  = tau->tauID("againstElectronVLooseMVA6");
    eventTuple().againstElectronVTightMVA6_1  = tau->tauID("againstElectronVTightMVA6");

    eventTuple().againstMuonLoose3_1          = tau->tauID("againstMuonLoose3");
    eventTuple().againstMuonTight3_1          = tau->tauID("againstMuonTight3");

    eventTuple().byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    eventTuple().byIsolationMVA3newDMwLTraw_1               = tau->tauID("byIsolationMVArun2v1DBnewDMwLTraw");
    eventTuple().byIsolationMVA3oldDMwLTraw_1               = tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
    eventTuple().byIsolationMVA3newDMwoLTraw_1              = default_value;
    eventTuple().byIsolationMVA3oldDMwoLTraw_1              = default_value;

    eventTuple().byVLooseIsolationMVArun2v1DBoldDMwLT_1     = tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byLooseIsolationMVArun2v1DBoldDMwLT_1      = tau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byMediumIsolationMVArun2v1DBoldDMwLT_1     = tau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byTightIsolationMVArun2v1DBoldDMwLT_1      = tau->tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byVTightIsolationMVArun2v1DBoldDMwLT_1     = tau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

    eventTuple().decayModeFindingOldDMs_1 = tau->tauID("decayModeFinding");

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau);
