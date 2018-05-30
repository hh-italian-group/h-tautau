/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muMu.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muMu::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuMu;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
    }

    const auto selectedMuons = CollectSignalMuons();
    cut(selectedMuons.size(), "muons");

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::MuMu::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::MuMu::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedMuons, selectedMuons, DeltaR_betweenSignalObjects, "H_mu_mu");
    cut(higgses.size(), "mu_mu_pair");

    cut(selected_higgs.GetFirstDaughter() < muonID::pfRelIso04 ||
        selected_higgs.GetSecondDaughter() < muonID::pfRelIso04, "iso_of_1_daughter");

    std::sort(higgses.begin(), higgses.end(), &HiggsComparitor<HiggsCandidate>);
    auto selected_higgs = higgses.front();
    if (selected_higgs.GetFirstDaughter().GetMomentum().Pt() < selected_higgs.GetSecondDaughter().GetMomentum().Pt())
        selected_higgs = HiggsCandidate(selected_higgs.GetSecondDaughter(), selected_higgs.GetFirstDaughter());

    if(applyTriggerMatch){
        triggerTools.SetTriggerMatchBits(selection.triggerResults, selected_higgs,
                                      cuts::H_tautau_2016::DeltaR_triggerMatch);
        cut(selection.triggerResults.AnyAcceptAndMatch(), "trigger_match");
    }

    selection.SetHiggsCandidate(selected_higgs);

    //Third-Lepton Veto
    const auto electronVetoCollection = CollectVetoElectrons();
    const auto muonVetoCollection = CollectVetoMuons({ &selection.higgs->GetFirstDaughter(),
                                                       &selection.higgs->GetSecondDaughter() });
    selection.electronVeto = electronVetoCollection.size();
    selection.muonVeto = muonVetoCollection.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);
    FillEventTuple(selection);

    if(eventEnergyScale == analysis::EventEnergyScale::Central)
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
}

std::vector<BaseTupleProducer::MuonCandidate> TupleProducer_muMu::CollectSignalMuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muMu::SelectSignalMuon, this, _1, _2);
    return CollectObjects("SignalMuons", base_selector, muons);
}

void TupleProducer_muMu::SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuMu::muonID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > cuts::hh_bbtautau_2016::MuMu::muonID::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);

    bool passMuonId = muon->isMediumMuon();
    if(productionMode == ProductionMode::hh)
        passMuonId = muon->isTightMuon(*primaryVertex);
    cut(passMuonId, "muonID");

//    if(productionMode == ProductionMode::hh)
//        cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());
}

void TupleProducer_muMu::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::MuMu);

    FillMuonLeg(1, selection.higgs->GetFirstDaughter());
    FillMuonLeg(2, selection.higgs->GetSecondDaughter());

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muMu);
