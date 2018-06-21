/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuTau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
    }

    // Signal-like leptons selection
    const auto selectedMuons = CollectSignalMuons();
    cut(selectedMuons.size(), "muons");
    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedMuons, selectedTaus, DeltaR_betweenSignalObjects, "H_mu_tau");
    cut(higgses.size(), "mu_tau_pair");

    std::sort(higgses.begin(), higgses.end(), &HiggsComparitor<HiggsCandidate>);
    auto selected_higgs = higgses.front();

    if(applyTriggerMatch){
        triggerTools.SetTriggerMatchBits(selection.triggerResults, selected_higgs,
                                      cuts::H_tautau_2016::DeltaR_triggerMatch);
        cut(selection.triggerResults.AnyAcceptAndMatch(), "trigger_match");
    }

    selection.SetHiggsCandidate(selected_higgs);

    //Third-Lepton Veto
    selection.other_electrons = CollectVetoElectrons();
    selection.other_muons = CollectVetoMuons({ &selection.higgs->GetFirstDaughter() });
    selection.electronVeto = selection.other_electrons.size();
    selection.muonVeto = selection.other_muons.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);
    FillEventTuple(selection);

    if(eventEnergyScale == analysis::EventEnergyScale::Central)
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
}

std::vector<BaseTupleProducer::MuonCandidate> TupleProducer_muTau::CollectSignalMuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalMuon, this, _1, _2);
    return CollectObjects("SignalMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_muTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_muTau::SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::muonID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    double pt_cut = pt;
    if(productionMode == ProductionMode::hh) {
        pt_cut = period == analysis::Period::Run2017 ? cuts::hh_bbtautau_2017::MuTau::muonID::pt : cuts::hh_bbtautau_2016::MuTau::muonID::pt;
    }
    else if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::MuTau::muonID::pt;
    else if (productionMode == ProductionMode::h_tt_sm) pt_cut = cuts::H_tautau_2016_sm::MuTau::muonID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    const double eta_cut  = productionMode == ProductionMode::h_tt_sm ? cuts::H_tautau_2016_sm::MuTau::muonID::eta : eta;
    cut(std::abs(p4.eta()) < eta_cut, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    if(productionMode == ProductionMode::hh){
        cut(muon->isTightMuon(*primaryVertex), "muonID");
        cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());
    }
    else cut(muon->isMediumMuon(), "muonID");

}

void TupleProducer_muTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    double pt_cut = pt;
    if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::MuTau::tauID::pt;
    cut(p4.Pt() > pt_cut, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const auto dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronVLooseMVA6") > againstElectronVLooseMVA6, "againstElectron");
        cut(tau->tauID("againstMuonTight3") > againstMuonTight3, "againstMuon");
        if(period == analysis::Period::Run2017) {
            cut(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5, "VVLooseIso");
        }
    }
}

void TupleProducer_muTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::MuTau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    const bool store_tauIds = !previous_selection || !selection.HaveSameSecondLegOrigin(*previous_selection);
    storageMode.SetPresence(EventPart::SecondTauIds, store_tauIds);
    eventTuple().storageMode = storageMode.Mode();

    FillMuonLeg(1, selection.higgs->GetFirstDaughter());
    FillTauLeg(2, selection.higgs->GetSecondDaughter(), store_tauIds);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muTau);
