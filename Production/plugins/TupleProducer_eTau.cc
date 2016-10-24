/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::ETau;

    SelectionResults selection;
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) cut(triggerTools.HaveTriggerFired(hltPaths), "trigger");

    //Di-Lepton Veto
    const auto z_electrons = CollectZelectrons();
    const auto z_electrons_candidates = FindCompatibleObjects(z_electrons, z_electrons, ZeeVeto::deltaR, "Z_e_e",0);
    selection.Zveto = z_electrons_candidates.size();

    // Signal-like leptons selection
    const auto selectedElectrons = CollectSignalElectrons();
    cut(selectedElectrons.size(), "electrons");
    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedElectrons, selectedTaus, DeltaR_betweenSignalObjects, "H_e_tau");
    cut(higgses.size(), "ele_tau_pair");

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
    using namespace cuts::H_tautau_2016::ETau::electronID;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    const double pt_cut = productionMode == ProductionMode::hh ? cuts::hh_bbtautau_2016::ETau::electronID::pt : pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_xy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_xy < dxy, "dxy", electron_xy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    const bool isTight = (*tight_id_decisions)[electron.getPtr()];
    cut(isTight, "electronMVATightID");
    if(productionMode != ProductionMode::hh) {
        const auto eleMissingHits =
                electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        cut(eleMissingHits <= missingHits, "missingHits", eleMissingHits);
        cut(electron->passConversionVeto(), "conversionVeto");
    } else {
        cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
    }
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::tauID;

    cut(true, "gt0_tau_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    cut(p4.Pt() > pt, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const auto dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronTightMVA6") > againstElectronTightMVA6, "againstElectron");
        cut(tau->tauID("againstMuonLoose3") > againstMuonLoose3, "againstMuon");
    }
}


void TupleProducer_eTau::FillEventTuple(const SelectionResults& selection)
{
    using namespace analysis;
    static constexpr float default_value = ntuple::DefaultFillValue<Float_t>();
    static constexpr int default_int_value = ntuple::DefaultFillValue<Int_t>();

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelID = static_cast<int>(analysis::Channel::ETau);

    // Leg 1, lepton
    const ElectronCandidate& electron = selection.higgs->GetFirstDaughter();
    eventTuple().p4_1     = LorentzVectorM(electron.GetMomentum());
    eventTuple().q_1      = electron.GetCharge();
    eventTuple().d0_1     = electron->gsfTrack()->dxy(primaryVertex->position());
    eventTuple().dZ_1     = electron->gsfTrack()->dz(primaryVertex->position());
    eventTuple().iso_1    = electron.GetIsolation();
    eventTuple().id_e_mva_nt_loose_1 = default_value;
    eventTuple().gen_match_1 = isMC ? gen_truth::genMatch(electron->p4(), *genParticles) : default_int_value;

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
