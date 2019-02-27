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

    // Signal-like leptons selection
    auto muons = CollectSignalMuons();
    cut(muons.size() > 1, "muons");

    std::sort(muons.begin(),muons.end(),&LeptonComparitor<MuonCandidate>);
    selection.muons.push_back(muons.at(0));
    selection.muons.push_back(muons.at(1));
    selection.other_muons = CollectVetoMuons({&muons.at(0),&muons.at(1)});
    selection.muonVeto = selection.other_muons.size();
    cut(!selection.muonVeto, "no_extra_muon");

    selection.other_electrons = CollectVetoElectrons();
    selection.electronVeto = selection.other_electrons.size();
    cut(!selection.electronVeto, "no_extra_ele");

    const double DeltaR_betweenSignalObjects = (productionMode == ProductionMode::hh ||
        productionMode == ProductionMode::tau_pog)
            ? cuts::hh_bbtautau_2016::MuMu::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::MuMu::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selection.muons, selection.muons, DeltaR_betweenSignalObjects, "H_mu_mu");
    cut(higgses.size(), "mu_mu_pair");

    for(size_t n = 0; n < higgses.size(); ++n){
        HiggsCandidate selected_higgs = higgses.at(n);
        if (higgses.at(n).GetFirstDaughter().GetMomentum().Pt() < higgses.at(n).GetSecondDaughter().GetMomentum().Pt())
            selected_higgs = HiggsCandidate(higgses.at(n).GetSecondDaughter(), higgses.at(n).GetFirstDaughter());

        cut(selected_higgs.GetFirstDaughter().GetIsolation() < muonID::pfRelIso04 ||
                selected_higgs.GetSecondDaughter().GetIsolation() < muonID::pfRelIso04, "iso_of_1_daughter");

        if(applyTriggerMatch){
            triggerTools.SetTriggerMatchBits(selection.triggerResults, selected_higgs,
                                          cuts::H_tautau_2016::DeltaR_triggerMatch);
        }
        selection.SetHiggsCandidate(selected_higgs);
        selection.higgses.push_back(selected_higgs);

        if(runSVfit)
            selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);

    }

    ApplyBaseSelection(selection);

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
}

void TupleProducer_muMu::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::MuMu);

    BaseTupleProducer::FillMuon(selection);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muMu);
